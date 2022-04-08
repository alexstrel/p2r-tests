#pragma once

#include <iterators.h>

namespace impl {
  /**
     Simple array object which mimics std::array
  */
  template <typename T, int n> struct array {
    using value_type = T;
    T data[n];

    constexpr T &operator[](int i) { return data[i]; }
    constexpr const T &operator[](int i) const { return data[i]; }
    constexpr int size() const { return n; }

    array() = default;
    array(const array<T, n> &) = default;
    array(array<T, n> &&) = default;

    array<T, n> &operator=(const array<T, n> &) = default;
    array<T, n> &operator=(array<T, n> &&) = default;
  };
} //impl


template <typename T, int N>
struct MPNX_ {
   impl::array<T,N> data;
   //basic accessors
   inline const T& operator[](const int idx) const {return data[idx];}
   inline T& operator[](const int idx) {return data[idx];}
};

using MP1I_    = MPNX_<int,   1 >;
using MP1F_    = MPNX_<float, 1 >;
using MP2F_    = MPNX_<float, 3 >;
using MP3F_    = MPNX_<float, 3 >;
using MP6F_    = MPNX_<float, 6 >;
using MP2x2SF_ = MPNX_<float, 3 >;
using MP3x3SF_ = MPNX_<float, 6 >;
using MP6x6SF_ = MPNX_<float, 21>;
using MP6x6F_  = MPNX_<float, 36>;
using MP3x3_   = MPNX_<float, 9 >;
using MP3x6_   = MPNX_<float, 18>;

template <typename T, typename Allocator, int n>
struct MPNX {
   using DataType = T;

   static constexpr int N    = n;

   const int nTrks;//note that bSize is a tuning parameter!
   const int nEvts;
   const int nLayers;

   std::vector<T, Allocator> data;

   MPNX() : nTrks(0), nEvts(0), nLayers(0), data(n){}

   MPNX(const int ntrks_, const int nevts_, const int nlayers_ = 1) :
      nTrks(ntrks_),
      nEvts(nevts_),
      nLayers(nlayers_),
      data(n*nTrks*nEvts*nLayers){
   }

   MPNX(const std::vector<T, Allocator> data_, const int ntrks_, const int nevts_, const int nlayers_ = 1) :
      nTrks(ntrks_),
      nEvts(nevts_),
      nLayers(nlayers_),
      data(data_) {
     if(data_.size() > n*nTrks*nEvts*nLayers) {std::cerr << "Incorrect dim parameters."; }
   }
};

using MP1I    = MPNX<int,  IntAllocator,   1 >;
using MP1F    = MPNX<float,FloatAllocator, 1 >;
using MP2F    = MPNX<float,FloatAllocator, 2 >;
using MP3F    = MPNX<float,FloatAllocator, 3 >;
using MP6F    = MPNX<float,FloatAllocator, 6 >;
using MP3x3   = MPNX<float,FloatAllocator, 9 >;
using MP3x6   = MPNX<float,FloatAllocator, 18>;
using MP2x2SF = MPNX<float,FloatAllocator, 3 >;
using MP3x3SF = MPNX<float,FloatAllocator, 6 >;
using MP6x6SF = MPNX<float,FloatAllocator, 21>;
using MP6x6F  = MPNX<float,FloatAllocator, 36>;

template <typename MPNTp, FieldOrder Order>
struct MPNXAccessor {
   typedef typename MPNTp::DataType T;

   static constexpr int n   = MPNTp::N;//matrix linear dim (total number of els)

   int nTrks;
   int nEvts;
   int nLayers;

   int NevtsNtrks;

   int stride;
   
   int thread_stride;

   T* data_; //accessor field only for the data access, not allocated here

   MPNXAccessor() = default;

   MPNXAccessor(const MPNTp &v) :
        nTrks(v.nTrks),
        nEvts(v.nEvts),
        nLayers(v.nLayers),
        NevtsNtrks(nEvts*nTrks),
        stride(Order == FieldOrder::P2R_TRACKBLK_EVENT_LAYER_MATIDX_ORDER ? nTrks*nEvts*nLayers  : nTrks*nEvts*n),
        thread_stride(Order == FieldOrder::P2R_TRACKBLK_EVENT_LAYER_MATIDX_ORDER ? stride  : NevtsNtrks),              
        data_(const_cast<T*>(v.data.data())){ }

   inline T& operator[](const int idx) const {return data_[idx];}

   inline T& operator()(const int mat_idx, const int trkev_idx, const int layer_idx) const {
     if      constexpr (Order == FieldOrder::P2R_TRACKBLK_EVENT_LAYER_MATIDX_ORDER)
       return data_[mat_idx*stride + layer_idx*NevtsNtrks + trkev_idx];//using defualt order batch id (the fastest) > track id > event id > layer id (the slowest)
     else //(Order == FieldOrder::P2R_TRACKBLK_EVENT_MATIDX_LAYER_ORDER)
       return data_[layer_idx*stride + mat_idx*NevtsNtrks + trkev_idx];
   }//i is the internal dof index

   inline T& operator()(const int thrd_idx, const int blk_offset) const { return data_[thrd_idx*thread_stride + blk_offset];}//

   inline int GetThreadOffset(const int thrd_idx, const int layer_idx = 0) const {
     if      constexpr (Order == FieldOrder::P2R_TRACKBLK_EVENT_LAYER_MATIDX_ORDER)
       return (layer_idx*NevtsNtrks + thrd_idx);//using defualt order batch id (the fastest) > track id > event id > layer id (the slowest)
     else //(Order == FieldOrder::P2R_TRACKBLK_EVENT_MATIDX_LAYER_ORDER)
       return (layer_idx*stride + thrd_idx);
   }
   
   inline void load(MPNX_<T, n>& dest, const int tid, const int layer = 0) const {
      auto tid_offset = GetThreadOffset(tid, layer);
#pragma unroll
      for(int id = 0; id < n; id++){
          dest[id] = this->operator()(id, tid_offset);
      }
      return;
   }
   inline void save(const MPNX_<T, n>& src, const int tid, const int layer = 0){
      auto tid_offset = GetThreadOffset(tid, layer); 
#pragma unroll
      for(int id = 0; id < n; id++){
        this->operator()(id, tid_offset) = src[id];
      }
      return;
   }  
  
};




