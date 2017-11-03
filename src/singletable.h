#ifndef CUCKOO_FILTER_SINGLE_TABLE_H_
#define CUCKOO_FILTER_SINGLE_TABLE_H_

#include <assert.h>

#include <sstream>

#include <xmmintrin.h>

#include "bitsutil.h"
#include "debug.h"
#include "printutil.h"

//#include "locklib.h"
#include "primitives.h"
#include "stats.h"
#include "cuckoofilter.h"



typedef qnode_t* mcs_lock_t;

typedef struct {
        mcs_lock_t    lock;
        uint32_t    pad1[128/4 - sizeof(mcs_lock_t)];
} spec_mcs_lock_t;

typedef struct {
        mcs_lock_t        mcs;
        uint8_t mode;
        char cv_lock_pad[128];
        pthread_mutex_t cv_lock;
} locklib_mutex_t;

#define MAX_RETRIES         10

static __thread qnode_t   my_node;
static __thread uint32_t  retries;

int locklib_mutex_lock(locklib_mutex_t *mutex, uint8_t mode)
{
        stats_begin();

        spec_mcs_lock_t *lock = (spec_mcs_lock_t *) mutex;
        uint32_t reason = 0;

        stats_count(locks);

speculative_path:
        XBEGIN(fallback_path, reason);
        if (lock->lock) XABORT(1);
        // if (lock->lock) {
        //   if (mutex->mode = 1) {
        //     XABORT(1);
        //   }
        // }
        return 0;

fallback_path:
        stats_count(contended);
        stats_count(reasons[reason & 63]);

        retries++;

        // lock taken?
        while (lock->lock)
                cpu_relax();

        if (retries < MAX_RETRIES)
                goto speculative_path;

        // Acquire lock in a standard manner
        my_node.locked = true;
        qnode_t *prev = __sync_lock_test_and_set(&lock->lock, &my_node);
        if (unlikely(prev != NULL)) {
                prev->next = &my_node;
                while (my_node.locked)
                        cpu_relax();
        }

        mutex->mode = mode;
        stats_count(failures);
        return 0;
}

int locklib_mutex_unlock(locklib_mutex_t *mutex)
{
        spec_mcs_lock_t *lock = (spec_mcs_lock_t *) mutex;

        if (XTEST()) {

                stats_count(speculative);
                XEND();
        } else {
                // Release lock in a standard manner
                qnode_t *last = my_node.next;
                if (last == NULL) {
                        if (likely(true == __sync_bool_compare_and_swap(&lock->lock, &my_node, NULL)))
                                return 0;

                        while ((last = my_node.next) == NULL)
                                cpu_relax();
                }

                my_node.next = NULL;
                last->locked = false;
        }

        retries = 0;
        return 0;
}

static inline uint8_t tas_uint8(volatile uint8_t *addr) {
    uint8_t oldval;
    __asm__ __volatile__("xchgb %0,%1"
            : "=q"(oldval), "=m"(*addr)
            : "0"((unsigned char) 0xff), "m"(*addr) : "memory");
    return (uint8_t) oldval;
}

#define TAS_U8(a) tas_uint8(a)

namespace cuckoofilter {

class SpinLock{
public:
  SpinLock(){
    spinlock_ = 0;
  }

  inline void Lock(){
    while(TAS_U8(&spinlock_)){
      _mm_pause();
    }
  }

  inline void Unlock(){
    __asm__ __volatile__("" ::: "memory");
    spinlock_ = 0;
  }

private:
  volatile uint8_t spinlock_;
};

enum ReturnCode {
  Succ = 0,
  Failed = 1,
  Exists = 2,
};

// the most naive table implementation: one huge bit array
template <size_t bits_per_tag>
class SingleTable {
  static const size_t kTagsPerBucket = 4;
  static const size_t kBytesPerBucket =
      (bits_per_tag * kTagsPerBucket + 7) >> 3;
  static const uint32_t kTagMask = (1ULL << bits_per_tag) - 1;

  struct Bucket {
    char bits_[kBytesPerBucket];
  } __attribute__((__packed__));

  struct SpinLock *locks;
  locklib_mutex_t *mutex;

  // using a pointer adds one more indirection
  Bucket *buckets_;
  size_t num_buckets_;

 public:
  explicit SingleTable(const size_t num) : num_buckets_(num) {
    buckets_ = new Bucket[num_buckets_];
    memset(buckets_, 0, kBytesPerBucket * num_buckets_);
    locks = new SpinLock[num_buckets_];
    mutex = new locklib_mutex_t[1];
    //locklib_mutex_init(mutex, NULL);
  }

  ~SingleTable() {
    delete[] buckets_;
    //locklib_mutex_destroy(mutex);
    delete[] mutex;
  }

  size_t NumBuckets() const {
    return num_buckets_;
  }

  size_t SizeInBytes() const {
    return kBytesPerBucket * num_buckets_;
  }

  size_t SizeInTags() const {
    return kTagsPerBucket * num_buckets_;
  }

  std::string Info() const {
    std::stringstream ss;
    ss << "SingleHashtable with tag size: " << bits_per_tag << " bits \n";
    ss << "\t\tAssociativity: " << kTagsPerBucket << "\n";
    ss << "\t\tTotal # of rows: " << num_buckets_ << "\n";
    ss << "\t\tTotal # slots: " << SizeInTags() << "\n";
    return ss.str();
  }

  // read tag from pos(i,j)
  inline uint32_t ReadTag(const size_t i, const size_t j) const {
    const char *p = buckets_[i].bits_;
    uint32_t tag;
    /* following code only works for little-endian */
    if (bits_per_tag == 2) {
      tag = *((uint8_t *)p) >> (j * 2);
    } else if (bits_per_tag == 4) {
      p += (j >> 1);
      tag = *((uint8_t *)p) >> ((j & 1) << 2);
    } else if (bits_per_tag == 8) {
      p += j;
      tag = *((uint8_t *)p);
    } else if (bits_per_tag == 12) {
      p += j + (j >> 1);
      tag = *((uint16_t *)p) >> ((j & 1) << 2);
    } else if (bits_per_tag == 16) {
      p += (j << 1);
      tag = *((uint16_t *)p);
    } else if (bits_per_tag == 32) {
      tag = ((uint32_t *)p)[j];
    }
    return tag & kTagMask;
  }

  // write tag to pos(i,j)
  inline void WriteTag(const size_t i, const size_t j, const uint32_t t) {
    char *p = buckets_[i].bits_;
    uint32_t tag = t & kTagMask;
    /* following code only works for little-endian */
    if (bits_per_tag == 2) {
      *((uint8_t *)p) |= tag << (2 * j);
    } else if (bits_per_tag == 4) {
      p += (j >> 1);
      if ((j & 1) == 0) {
        *((uint8_t *)p) &= 0xf0;
        *((uint8_t *)p) |= tag;
      } else {
        *((uint8_t *)p) &= 0x0f;
        *((uint8_t *)p) |= (tag << 4);
      }
    } else if (bits_per_tag == 8) {
      ((uint8_t *)p)[j] = tag;
    } else if (bits_per_tag == 12) {
      p += (j + (j >> 1));
      if ((j & 1) == 0) {
        ((uint16_t *)p)[0] &= 0xf000;
        ((uint16_t *)p)[0] |= tag;
      } else {
        ((uint16_t *)p)[0] &= 0x000f;
        ((uint16_t *)p)[0] |= (tag << 4);
      }
    } else if (bits_per_tag == 16) {
      ((uint16_t *)p)[j] = tag;
    } else if (bits_per_tag == 32) {
      ((uint32_t *)p)[j] = tag;
    }
  }

  inline bool FindTagInBuckets(const size_t i1, const size_t i2,
                               const uint32_t tag) const {
    const char *p1 = buckets_[i1].bits_;
    const char *p2 = buckets_[i2].bits_;

    uint64_t v1 = *((uint64_t *)p1);
    uint64_t v2 = *((uint64_t *)p2);

    // caution: unaligned access & assuming little endian
    if (bits_per_tag == 4 && kTagsPerBucket == 4) {
      return hasvalue4(v1, tag) || hasvalue4(v2, tag);
    } else if (bits_per_tag == 8 && kTagsPerBucket == 4) {
      return hasvalue8(v1, tag) || hasvalue8(v2, tag);
    } else if (bits_per_tag == 12 && kTagsPerBucket == 4) {
      return hasvalue12(v1, tag) || hasvalue12(v2, tag);
    } else if (bits_per_tag == 16 && kTagsPerBucket == 4) {
      return hasvalue16(v1, tag) || hasvalue16(v2, tag);
    } else {
      for (size_t j = 0; j < kTagsPerBucket; j++) {
        if ((ReadTag(i1, j) == tag) || (ReadTag(i2, j) == tag)) {
          return true;
        }
      }
      return false;
    }
  }

#ifdef HTM_LOCK 
  #define FIND_LOCK     \
          locklib_mutex_lock(mutex, 0)
  #define FIND_UNLOCK   \
          locklib_mutex_unlock(mutex)
  #define UPDATE_LOCK   \
          locklib_mutex_lock(mutex, 1)
  #define UPDATE_UNLOCK \
          locklib_mutex_unlock(mutex)
#else
  #define FIND_LOCK           \
          locks[i1].Lock();   \
          locks[i2].Lock();
  #define FIND_UNLOCK         \
          locks[i2].Unlock(); \
          locks[i1].Unlock();
  #define UPDATE_LOCK         \
          locks[i].Lock();
  #define UPDATE_UNLOCK       \
          locks[i].Unlock();
#endif  

  inline bool ConFindTagInBuckets(const size_t i1, const size_t i2,
                               const uint32_t tag) const {
    FIND_LOCK;
    const char *p1 = buckets_[i1].bits_;
    const char *p2 = buckets_[i2].bits_;

    uint64_t v1 = *((uint64_t *)p1);
    uint64_t v2 = *((uint64_t *)p2);
    
    // caution: unaligned access & assuming little endian
    if (bits_per_tag == 4 && kTagsPerBucket == 4) {
      bool ret = hasvalue4(v1, tag) || hasvalue4(v2, tag);
      FIND_UNLOCK;
      return ret;
    } else if (bits_per_tag == 8 && kTagsPerBucket == 4) {
      bool ret = hasvalue8(v1, tag) || hasvalue8(v2, tag);
      FIND_UNLOCK;
      return ret;
    } else if (bits_per_tag == 12 && kTagsPerBucket == 4) {
      // uint64_t v1old = v1;
      // uint64_t v2old = v2;
      // bool ret1 = hasvalue12(v1, tag);
      // bool ret2 = hasvalue12(v2, tag);
      // FIND_LOCK;
      // bool ret;
      // if (v1 == v1old || v2 == v2old) {
      //   ret = ret1 || ret2;
      // } else {
      //   ret = hasvalue12(v1, tag) || hasvalue12(v2, tag);
      // }
      bool ret = hasvalue12(v1, tag) || hasvalue12(v2, tag);
      FIND_UNLOCK;
      return ret;
    } else if (bits_per_tag == 16 && kTagsPerBucket == 4) {
      bool ret = hasvalue16(v1, tag) || hasvalue16(v2, tag);
      FIND_UNLOCK;
      return ret;
    } else {
      for (size_t j = 0; j < kTagsPerBucket; j++) {
        if ((ReadTag(i1, j) == tag) || (ReadTag(i2, j) == tag)) {
          FIND_UNLOCK;
          return true;
        }
      }
      FIND_UNLOCK;
      return false;
    }
  }

  inline bool FindTagInBucket(const size_t i, const uint32_t tag) const {
    // caution: unaligned access & assuming little endian
    if (bits_per_tag == 4 && kTagsPerBucket == 4) {
      const char *p = buckets_[i].bits_;
      uint64_t v = *(uint64_t *)p;  // uint16_t may suffice
      return hasvalue4(v, tag);
    } else if (bits_per_tag == 8 && kTagsPerBucket == 4) {
      const char *p = buckets_[i].bits_;
      uint64_t v = *(uint64_t *)p;  // uint32_t may suffice
      return hasvalue8(v, tag);
    } else if (bits_per_tag == 12 && kTagsPerBucket == 4) {
      const char *p = buckets_[i].bits_;
      uint64_t v = *(uint64_t *)p;
      return hasvalue12(v, tag);
    } else if (bits_per_tag == 16 && kTagsPerBucket == 4) {
      const char *p = buckets_[i].bits_;
      uint64_t v = *(uint64_t *)p;
      return hasvalue16(v, tag);
    } else {
      for (size_t j = 0; j < kTagsPerBucket; j++) {
        if (ReadTag(i, j) == tag) {
          return true;
        }
      }
      return false;
    }
  }

  inline bool ConDeleteTagFromBucket1(const size_t i, const uint32_t tag) {
    const char *oldp = buckets_[i].bits_;
    
    //uint32_t tags[4];
    int k = -1;
    for (size_t j = 0; j < kTagsPerBucket; j++) {
      //tags[j] = ReadTag(i, j);
      if (ReadTag(i, j) == tag) {
        k = j;
        break;
      }
    }
    if (k == -1)
      return false;
    oldp += k + (k >> 1);
    UPDATE_LOCK;
    const char *newp = buckets_[i].bits_;
    newp += k + (k >> 1);
    if (*((uint16_t*)oldp) == *((uint16_t*)newp)) {
      WriteTag(i, k, 0);
      UPDATE_UNLOCK;
      return true;
    
    } else {
    for (size_t j = 0; j < kTagsPerBucket; j++) {
      if (ReadTag(i, j) == tag) {
        assert(FindTagInBucket(i, tag) == true);
        WriteTag(i, j, 0);
        UPDATE_UNLOCK;
        return true;
      }
    }
  }
    UPDATE_UNLOCK;
    return false;
  }

  inline bool ConDeleteTagFromBucket(const size_t i, const uint32_t tag) {
   
    UPDATE_LOCK;    
    for (size_t j = 0; j < kTagsPerBucket; j++) {
      if (ReadTag(i, j) == tag) {
        assert(FindTagInBucket(i, tag) == true);
        WriteTag(i, j, 0);
        UPDATE_UNLOCK;
        return true;
      }
    }
    UPDATE_UNLOCK;
    return false;
  }

  inline bool DeleteTagFromBucket(const size_t i, const uint32_t tag) {
    for (size_t j = 0; j < kTagsPerBucket; j++) {
      if (ReadTag(i, j) == tag) {
        assert(FindTagInBucket(i, tag) == true);
        WriteTag(i, j, 0);
        return true;
      }
    }
    return false;
  }

  inline ReturnCode ConInsertTagToBucket(const size_t i, const uint32_t tag,
                                const bool kickout, uint32_t &oldtag) {
    UPDATE_LOCK;
    if (FindTagInBucket(i, tag) == true) {
      UPDATE_UNLOCK;
      return Exists;
    }

    for (size_t j = 0; j < kTagsPerBucket; j++) {
      if (ReadTag(i, j) == 0) {
        WriteTag(i, j, tag);
        UPDATE_UNLOCK;
        return Succ;
      }
    }
    UPDATE_UNLOCK;

    if (kickout) {
      size_t r = rand() % kTagsPerBucket;
      UPDATE_LOCK;
      oldtag = ReadTag(i, r);
      WriteTag(i, r, tag);
      UPDATE_UNLOCK;
    }
    //printf("insert failed.\n");
    return Failed;
  }

  inline bool InsertTagToBucket(const size_t i, const uint32_t tag,
                                const bool kickout, uint32_t &oldtag) {
    for (size_t j = 0; j < kTagsPerBucket; j++) {
      if (ReadTag(i, j) == 0) {
        WriteTag(i, j, tag);
        return true;
      }
    }
    if (kickout) {
      size_t r = rand() % kTagsPerBucket;
      oldtag = ReadTag(i, r);
      WriteTag(i, r, tag);
    }
    return false;
  }

  inline size_t NumTagsInBucket(const size_t i) const {
    size_t num = 0;
    for (size_t j = 0; j < kTagsPerBucket; j++) {
      if (ReadTag(i, j) != 0) {
        num++;
      }
    }
    return num;
  }

  inline size_t TotalItems() const {
    size_t num = 0;
    for (size_t i = 0; i < num_buckets_; i++) {
      for (size_t j = 0; j < kTagsPerBucket; j++)
        if (ReadTag(i, j) != 0)
          num++;
    }
    return num;
  }
};
}  // namespace cuckoofilter
#endif  // CUCKOO_FILTER_SINGLE_TABLE_H_
