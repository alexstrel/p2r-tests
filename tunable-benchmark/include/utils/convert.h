#pragma once

#include <algorithm>
#include <vector>
#include <memory>
#include <numeric>
#include <execution>

#include <p2r_enum.h>

//Forward declarations
struct MPTRK_;
struct MPHIT_;
struct MPTRK;
struct MPHIT;

template<FieldOrder order, typename MPTRKAllocator, ConversionType convers_tp>
void convertTracks(std::vector<MPTRK_, MPTRKAllocator> &external_order_data, MPTRK* internal_order_data) {
  //create an accessor field:
  std::unique_ptr<MPTRKAccessor<order>> ind(new MPTRKAccessor<order>(*internal_order_data));
  // store in element order for bunches of bsize matrices (a la matriplex)
  const int outer_loop_range = nevts*ntrks;
  //
  auto policy = std::execution::par_unseq;
  //
  std::for_each(policy, 
		impl::counting_iterator(0),
                impl::counting_iterator(outer_loop_range),
                [=, exd_ = external_order_data.data(), &ind_ = *ind] (const auto tid) {
                  {
                  //const int l = it+ib*bsize+ie*ntrks*bsize;
                    //par
    	            for (int ip=0;ip<6;++ip) {
    	              if constexpr (convers_tp == ConversionType::P2R_CONVERT_FROM_INTERNAL_ORDER)
    	                exd_[tid].par.data[ip] = ind_.par(ip, tid, 0);
    	              else
    	                ind_.par(ip, tid, 0) = exd_[tid].par.data[ip];  
    	            }
    	            //cov
    	            for (int ip=0;ip<21;++ip) {
    	              if constexpr (convers_tp == ConversionType::P2R_CONVERT_FROM_INTERNAL_ORDER)
    	                exd_[tid].cov.data[ip] = ind_.cov(ip, tid, 0);
    	              else
    	                ind_.cov(ip, tid, 0) = exd_[tid].cov.data[ip];
    	            }
    	            //q
    	            if constexpr (convers_tp == ConversionType::P2R_CONVERT_FROM_INTERNAL_ORDER)
    	              exd_[tid].q.data[0] = ind_.q(0, tid, 0);//fixme check
    	            else
    	              ind_.q(0, tid, 0) = exd_[tid].q.data[0];
                  }
                });
   //
   return;
}


template<FieldOrder order, typename MPHITAllocator, ConversionType convers_tp>
void convertHits(std::vector<MPHIT_, MPHITAllocator> &external_order_data, MPHIT* internal_oder_data) {
  //create an accessor field:
  std::unique_ptr<MPHITAccessor<order>> ind(new MPHITAccessor<order>(*internal_oder_data));
  // store in element order for bunches of bsize matrices (a la matriplex)
  const int outer_loop_range = nevts*ntrks;
  //
  auto policy = std::execution::par_unseq;
  //
  std::for_each(policy,
		impl::counting_iterator(0),
                impl::counting_iterator(outer_loop_range),
                [=, exd_ = external_order_data.data(), &ind_ = *ind] (const auto tid) {
                   //  
                   for(int layer=0; layer<nlayer; ++layer) {  
                     {
                       //pos
                       for (int ip=0;ip<3;++ip) {
                         if constexpr (convers_tp == ConversionType::P2R_CONVERT_FROM_INTERNAL_ORDER)
                           exd_[layer+nlayer*tid].pos.data[ip] = ind_.pos(ip, tid, layer);
                         else
                           ind_.pos(ip, tid, layer) = exd_[layer+nlayer*tid].pos.data[ip];
                       }
                       //cov
                       for (int ip=0;ip<6;++ip) {
                         if constexpr (convers_tp == ConversionType::P2R_CONVERT_FROM_INTERNAL_ORDER)
                           exd_[layer+nlayer*tid].cov.data[ip] = ind_.cov(ip, tid, layer);
                         else
                           ind_.cov(ip, tid, layer) = exd_[layer+nlayer*tid].cov.data[ip];
                       }
                     } 
                  }
               });
  
  return;
}




