/*
 nvc++ -cuda -O2 -std=c++20 -stdpar=gpu -gpu=cc86 -gpu=managed -gpu=fma -gpu=fastmath -gpu=autocollapse -gpu=loadcache:L1 -gpu=unroll ./src/propagate-tor-test_cuda_hybrid.cpp  -o ./"propagate_nvcpp_cuda"
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <sys/time.h>

#include <cassert>

#include <algorithm>
#include <vector>
#include <memory>
#include <numeric>
#include <execution>
#include <random>

//#ifdef _NVHPC_CUDA
#include <nv/target>
//#endif

#include <iterators.h>

const std::array<int, 36> SymOffsets66{0, 1, 3, 6, 10, 15, 1, 2, 4, 7, 11, 16, 3, 4, 5, 8, 12, 17, 6, 7, 8, 9, 13, 18, 10, 11, 12, 13, 14, 19, 15, 16, 17, 18, 19, 20};

struct ATRK {
  std::array<float,6> par;
  std::array<float,21> cov;
  int q;
};

struct AHIT {
  std::array<float,3> pos;
  std::array<float,6> cov;
};

constexpr int iparX     = 0;
constexpr int iparY     = 1;
constexpr int iparZ     = 2;
constexpr int iparIpt   = 3;
constexpr int iparPhi   = 4;
constexpr int iparTheta = 5;


struct MPTRK_ {
  MP6F_    par;
  MP6x6SF_ cov;
  MP1I_    q;

  //  MP22I   hitidx;
};

struct MPHIT_ {
  MP3F_    pos;
  MP3x3SF_ cov;
};

using MPTRKAllocator = impl::UVMAllocator<MPTRK_>;
using MPHITAllocator = impl::UVMAllocator<MPHIT_>;

struct MPTRK {
  MP6F    par;
  MP6x6SF cov;
  MP1I    q;

  MPTRK() : par(), cov(), q() {}
  MPTRK(const int ntrks_, const int nevts_) : par(ntrks_, nevts_), cov(ntrks_, nevts_), q(ntrks_, nevts_) {}
  //  MP22I   hitidx;
};

template <FieldOrder Order>
struct MPTRKAccessor {
  using MP6FAccessor   = MPNXAccessor<MP6F,    Order>;
  using MP6x6SFAccessor= MPNXAccessor<MP6x6SF, Order>;
  using MP1IAccessor   = MPNXAccessor<MP1I,    Order>;

  MP6FAccessor    par;
  MP6x6SFAccessor cov;
  MP1IAccessor    q;

  MPTRKAccessor() : par(), cov(), q() {}
  MPTRKAccessor(const MPTRK &in) : par(in.par), cov(in.cov), q(in.q) {}
  
  inline void load(MPTRK_ &dst, const int tid, const int layer = 0) const {
    this->par.load(dst.par, tid, layer);
    this->cov.load(dst.cov, tid, layer);
    this->q.load(dst.q, tid, layer);
    
    return;
  }
  
  inline void save(MPTRK_ &src, const int tid, const int layer = 0) {
    this->par.save(src.par, tid, layer);
    this->cov.save(src.cov, tid, layer);
    this->q.save(src.q, tid, layer);
    
    return;
  }
};

struct MPHIT {
  MP3F    pos;
  MP3x3SF cov;

  MPHIT() : pos(), cov(){}
  MPHIT(const int ntrks_, const int nevts_, const int nlayers_) : pos(ntrks_, nevts_, nlayers_), cov(ntrks_, nevts_, nlayers_) {}
};

template <FieldOrder Order>
struct MPHITAccessor {
  using MP3FAccessor   = MPNXAccessor<MP3F,    Order>;
  using MP3x3SFAccessor= MPNXAccessor<MP3x3SF, Order>;

  MP3FAccessor    pos;
  MP3x3SFAccessor cov;

  MPHITAccessor() : pos(), cov() {}
  MPHITAccessor(const MPHIT &in) : pos(in.pos), cov(in.cov) {}
  
  void load(MPHIT_ &dst, const int tid, const int layer = 0) const {
    this->pos.load(dst.pos, tid, layer);
    this->cov.load(dst.cov, tid, layer);
    
    return;
  }
  
  void save(MPHIT_ &src, const int tid, const int layer = 0) {
    this->pos.save(src.pos, tid, layer);
    this->cov.save(src.cov, tid, layer);
    
    return;
  } 
};


///////////////////////////////////////
//Gen. utils

float randn(float mu, float sigma) {
  float U1, U2, W, mult;
  static float X1, X2;
  static int call = 0;
  if (call == 1) {
    call = !call;
    return (mu + sigma * (float) X2);
  } do {
    U1 = -1 + ((float) rand () / RAND_MAX) * 2;
    U2 = -1 + ((float) rand () / RAND_MAX) * 2;
    W = pow (U1, 2) + pow (U2, 2);
  }
  while (W >= 1 || W == 0); 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult; 
  call = !call; 
  return (mu + sigma * (float) X1);
}


template<typename MPTRKAllocator>
void prepareTracks(std::vector<MPTRK_, MPTRKAllocator> &trcks, ATRK &inputtrk) {
  //
  for (int ie=0;ie<nevts;++ie) {
    for (int ib=0;ib<ntrks;++ib) {
      {
	      //par
	      for (int ip=0;ip<6;++ip) {
	        trcks[ib + ntrks*ie].par.data[ip] = (1+smear*randn(0,1))*inputtrk.par[ip];
	      }
	      //cov, scale by factor 100
	      for (int ip=0;ip<21;++ip) {
	        trcks[ib + ntrks*ie].cov.data[ip] = (1+smear*randn(0,1))*inputtrk.cov[ip]*100;
	      }
	      //q
	      trcks[ib + ntrks*ie].q.data[0] = inputtrk.q;//can't really smear this or fit will be wrong
      }
    }
  }
  //
  return;
}

template<typename MPHITAllocator>
void prepareHits(std::vector<MPHIT_, MPHITAllocator> &hits, std::vector<AHIT>& inputhits) {
  // store in element order for bunches of bsize matrices (a la matriplex)
  for (int lay=0;lay<nlayer;++lay) {

    int mylay = lay;
    if (lay>=inputhits.size()) {
      // int wraplay = inputhits.size()/lay;
      exit(1);
    }
    AHIT& inputhit = inputhits[mylay];

    for (int ie=0;ie<nevts;++ie) {
      for (int ib=0;ib<ntrks;++ib) {
        {
        	//pos
        	for (int ip=0;ip<3;++ip) {
        	  hits[lay+nlayer*(ib + ntrks*ie)].pos.data[ip] = (1+smear*randn(0,1))*inputhit.pos[ip];
        	}
        	//cov
        	for (int ip=0;ip<6;++ip) {
        	  hits[lay+nlayer*(ib + ntrks*ie)].cov.data[ip] = (1+smear*randn(0,1))*inputhit.cov[ip];
        	}
        }
      }
    }
  }
  return;
}


//////////////////////////////////////////////////////////////////////////////////////
// Aux utils 
MPTRK_* bTk(MPTRK_* tracks, int ev, int ib) {
  return &(tracks[ib + ntrks*ev]);
}

const MPTRK_* bTk(const MPTRK_* tracks, int ev, int ib) {
  return &(tracks[ib + ntrks*ev]);
}

float q(const MP1I_* bq, int it){
  return (*bq).data[0];
}
//
float par(const MP6F_* bpars, int it, int ipar){
  return (*bpars).data[it + ipar];
}
float x    (const MP6F_* bpars, int it){ return par(bpars, it, 0); }
float y    (const MP6F_* bpars, int it){ return par(bpars, it, 1); }
float z    (const MP6F_* bpars, int it){ return par(bpars, it, 2); }
float ipt  (const MP6F_* bpars, int it){ return par(bpars, it, 3); }
float phi  (const MP6F_* bpars, int it){ return par(bpars, it, 4); }
float theta(const MP6F_* bpars, int it){ return par(bpars, it, 5); }
//
float par(const MPTRK_* btracks, int it, int ipar){
  return par(&(*btracks).par,it,ipar);
}
float x    (const MPTRK_* btracks, int it){ return par(btracks, it, 0); }
float y    (const MPTRK_* btracks, int it){ return par(btracks, it, 1); }
float z    (const MPTRK_* btracks, int it){ return par(btracks, it, 2); }
float ipt  (const MPTRK_* btracks, int it){ return par(btracks, it, 3); }
float phi  (const MPTRK_* btracks, int it){ return par(btracks, it, 4); }
float theta(const MPTRK_* btracks, int it){ return par(btracks, it, 5); }
//
float par(const MPTRK_* tracks, int ev, int tk, int ipar){
  int ib = tk;
  const MPTRK_* btracks = bTk(tracks, ev, ib);
  int it = 0;
  return par(btracks, it, ipar);
}
float x    (const MPTRK_* tracks, int ev, int tk){ return par(tracks, ev, tk, 0); }
float y    (const MPTRK_* tracks, int ev, int tk){ return par(tracks, ev, tk, 1); }
float z    (const MPTRK_* tracks, int ev, int tk){ return par(tracks, ev, tk, 2); }
float ipt  (const MPTRK_* tracks, int ev, int tk){ return par(tracks, ev, tk, 3); }
float phi  (const MPTRK_* tracks, int ev, int tk){ return par(tracks, ev, tk, 4); }
float theta(const MPTRK_* tracks, int ev, int tk){ return par(tracks, ev, tk, 5); }
//

const MPHIT_* bHit(const MPHIT_* hits, int ev, int ib) {
  return &(hits[ib + ntrks*ev]);
}
const MPHIT_* bHit(const MPHIT_* hits, int ev, int ib,int lay) {
return &(hits[lay + (ib*nlayer) +(ev*nlayer*ntrks)]);
}
//
float Pos(const MP3F_* hpos, int it, int ipar){
  return (*hpos).data[it + ipar];
}
float x(const MP3F_* hpos, int it)    { return Pos(hpos, it, 0); }
float y(const MP3F_* hpos, int it)    { return Pos(hpos, it, 1); }
float z(const MP3F_* hpos, int it)    { return Pos(hpos, it, 2); }
//
float Pos(const MPHIT_* hits, int it, int ipar){
  return Pos(&(*hits).pos,it,ipar);
}
float x(const MPHIT_* hits, int it)    { return Pos(hits, it, 0); }
float y(const MPHIT_* hits, int it)    { return Pos(hits, it, 1); }
float z(const MPHIT_* hits, int it)    { return Pos(hits, it, 2); }
//
float Pos(const MPHIT_* hits, int ev, int tk, int ipar){
  int ib = tk;
  const MPHIT_* bhits = bHit(hits, ev, ib);
  int it = 0;
  return Pos(bhits,it,ipar);
}
float x(const MPHIT_* hits, int ev, int tk)    { return Pos(hits, ev, tk, 0); }
float y(const MPHIT_* hits, int ev, int tk)    { return Pos(hits, ev, tk, 1); }
float z(const MPHIT_* hits, int ev, int tk)    { return Pos(hits, ev, tk, 2); }


////////////////////////////////////////////////////////////////////////
///MAIN compute kernels


template <FieldOrder order = FieldOrder::P2R_TRACKBLK_EVENT_LAYER_MATIDX_ORDER, bool grid_stride = true>
__global__ void launch_p2r_kernels(MPTRKAccessor<order> &obtracksAcc, MPTRKAccessor<order> &btracksAcc, MPHITAccessor<order> &bhitsAcc, const int length){
   auto i = threadIdx.x + blockIdx.x * blockDim.x;

   MPTRK_ btracks;
   MPTRK_ obtracks;
   MPHIT_ bhits;

   while (i < length) {
     //
     btracksAcc.load(btracks, i);
     for(int layer=0; layer<nlayer; ++layer) {  
       //
       bhitsAcc.load(bhits, i, layer);
       //
       propagateToR(btracks.cov, btracks.par, btracks.q, bhits.pos, obtracks.cov, obtracks.par);
       KalmanUpdate(obtracks.cov, obtracks.par, bhits.cov, bhits.pos);
       //
     }
     //
     obtracksAcc.save(obtracks, i);
     
     if (grid_stride)
       i += gridDim.x * blockDim.x;
     else
       break;
  }
  return;
}


int main (int argc, char* argv[]) {

   #include "input_track.h"

   std::vector<AHIT> inputhits{inputhit21,inputhit20,inputhit19,inputhit18,inputhit17,inputhit16,inputhit15,inputhit14,
                               inputhit13,inputhit12,inputhit11,inputhit10,inputhit09,inputhit08,inputhit07,inputhit06,
                               inputhit05,inputhit04,inputhit03,inputhit02,inputhit01,inputhit00};

   printf("track in pos: x=%f, y=%f, z=%f, r=%f, pt=%f, phi=%f, theta=%f \n", inputtrk.par[0], inputtrk.par[1], inputtrk.par[2],
	  sqrtf(inputtrk.par[0]*inputtrk.par[0] + inputtrk.par[1]*inputtrk.par[1]),
	  1./inputtrk.par[3], inputtrk.par[4], inputtrk.par[5]);
   printf("track in cov: xx=%.2e, yy=%.2e, zz=%.2e \n", inputtrk.cov[SymOffsets66[0]],
	                                       inputtrk.cov[SymOffsets66[(1*6+1)]],
	                                       inputtrk.cov[SymOffsets66[(2*6+2)]]);
   for (int lay=0; lay<nlayer; lay++){
     printf("hit in layer=%lu, pos: x=%f, y=%f, z=%f, r=%f \n", lay, inputhits[lay].pos[0], inputhits[lay].pos[1], inputhits[lay].pos[2], sqrtf(inputhits[lay].pos[0]*inputhits[lay].pos[0] + inputhits[lay].pos[1]*inputhits[lay].pos[1]));
   }
   
   printf("produce nevts=%i ntrks=%i smearing by=%f \n", nevts, ntrks, smear);
   printf("NITER=%d\n", NITER);

   long setup_start, setup_stop;
   struct timeval timecheck;

   constexpr auto order = FieldOrder::P2R_TRACKBLK_EVENT_LAYER_MATIDX_ORDER;

   using MPTRKAccessorTp = MPTRKAccessor<order>;
   using MPHITAccessorTp = MPHITAccessor<order>;

   impl::UVMAllocator<MPTRKAccessorTp> mptrk_uvm_alloc;
   impl::UVMAllocator<MPHITAccessorTp> mphit_uvm_alloc;

   gettimeofday(&timecheck, NULL);
   setup_start = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

   std::unique_ptr<MPTRK> trcksPtr(new MPTRK(ntrks, nevts));
   auto trcksAccPtr = std::allocate_shared<MPTRKAccessorTp>(mptrk_uvm_alloc, *trcksPtr);
   //
   std::unique_ptr<MPHIT> hitsPtr(new MPHIT(ntrks, nevts, nlayer));
   auto hitsAccPtr = std::allocate_shared<MPHITAccessorTp>(mphit_uvm_alloc, *hitsPtr);
   //
   std::unique_ptr<MPTRK> outtrcksPtr(new MPTRK(ntrks, nevts));
   auto outtrcksAccPtr = std::allocate_shared<MPTRKAccessorTp>(mptrk_uvm_alloc, *outtrcksPtr);
   //
   using hostmptrk_allocator = std::allocator<MPTRK_>;
   using hostmphit_allocator = std::allocator<MPHIT_>;

   std::vector<MPTRK_, hostmptrk_allocator > trcks(nevts*ntrks); 
   prepareTracks<hostmptrk_allocator>(trcks, inputtrk);
   //
   std::vector<MPHIT_, hostmphit_allocator> hits(nlayer*nevts*ntrks);
   prepareHits<hostmphit_allocator>(hits, inputhits);
   //
   std::vector<MPTRK_, hostmptrk_allocator> outtrcks(nevts*ntrks);
   
   convertHits<order, hostmphit_allocator, ConversionType::P2R_CONVERT_TO_INTERNAL_ORDER>(hits,     hitsPtr.get());
   convertTracks<order, hostmptrk_allocator, ConversionType::P2R_CONVERT_TO_INTERNAL_ORDER>(trcks,    trcksPtr.get());
   convertTracks<order, hostmptrk_allocator, ConversionType::P2R_CONVERT_TO_INTERNAL_ORDER>(outtrcks, outtrcksPtr.get());

   cudaDeviceSynchronize();

   gettimeofday(&timecheck, NULL);
   setup_stop = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

   printf("done preparing!\n");

   printf("Size of struct MPTRK trk[] = %ld\n", nevts*ntrks*sizeof(MPTRK));
   printf("Size of struct MPTRK outtrk[] = %ld\n", nevts*ntrks*sizeof(MPTRK));
   printf("Size of struct struct MPHIT hit[] = %ld\n", nevts*ntrks*sizeof(MPHIT));

   const int phys_length      = nevts*ntrks;
   const int outer_loop_range = phys_length;
   //
   auto wall_start = std::chrono::high_resolution_clock::now();

   dim3 blocks(threadsperblock, 1, 1);
   dim3 grid(((outer_loop_range + threadsperblock - 1)/ threadsperblock),1,1);

   for(int itr=0; itr<NITER; itr++) {

     launch_p2r_kernels<<<grid, blocks>>>(*outtrcksAccPtr, *trcksAccPtr, *hitsAccPtr, phys_length);

   } //end of itr loop

   cudaDeviceSynchronize();

   auto wall_stop = std::chrono::high_resolution_clock::now();

   auto wall_diff = wall_stop - wall_start;
   auto wall_time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(wall_diff).count()) / 1e6;   

   printf("setup time time=%f (s)\n", (setup_stop-setup_start)*0.001);
   printf("done ntracks=%i tot time=%f (s) time/trk=%e (s)\n", nevts*ntrks*int(NITER), wall_time, wall_time/(nevts*ntrks*int(NITER)));
   printf("formatted %i %i %i %i %i %f 0 %f %i\n",int(NITER),nevts, ntrks, 1, ntrks, wall_time, (setup_stop-setup_start)*0.001, -1);

   convertTracks<order, hostmptrk_allocator, ConversionType::P2R_CONVERT_FROM_INTERNAL_ORDER>(outtrcks, outtrcksPtr.get());
   auto outtrk = outtrcks.data();

   int nnans = 0, nfail = 0;
   float avgx = 0, avgy = 0, avgz = 0, avgr = 0;
   float avgpt = 0, avgphi = 0, avgtheta = 0;
   float avgdx = 0, avgdy = 0, avgdz = 0, avgdr = 0;

   for (int ie=0;ie<nevts;++ie) {
     for (int it=0;it<ntrks;++it) {
       float x_ = x(outtrk,ie,it);
       float y_ = y(outtrk,ie,it);
       float z_ = z(outtrk,ie,it);
       float r_ = sqrtf(x_*x_ + y_*y_);
       float pt_ = std::abs(1./ipt(outtrk,ie,it));
       float phi_ = phi(outtrk,ie,it);
       float theta_ = theta(outtrk,ie,it);
       float hx_ = inputhits[nlayer-1].pos[0];
       float hy_ = inputhits[nlayer-1].pos[1];
       float hz_ = inputhits[nlayer-1].pos[2];
       float hr_ = sqrtf(hx_*hx_ + hy_*hy_);
 
       if (std::isfinite(x_)==false ||
          std::isfinite(y_)==false ||
          std::isfinite(z_)==false ||
          std::isfinite(pt_)==false ||
          std::isfinite(phi_)==false ||
          std::isfinite(theta_)==false
          ) {
        nnans++;
        continue;
       }
       if (fabs( (x_-hx_)/hx_ )>1. ||
	   fabs( (y_-hy_)/hy_ )>1. ||
	   fabs( (z_-hz_)/hz_ )>1.) {
	 nfail++;
	 continue;
       }

       avgpt += pt_;
       avgphi += phi_;
       avgtheta += theta_;
       avgx += x_;
       avgy += y_;
       avgz += z_;
       avgr += r_;
       avgdx += (x_-hx_)/x_;
       avgdy += (y_-hy_)/y_;
       avgdz += (z_-hz_)/z_;
       avgdr += (r_-hr_)/r_;
       //if((it+ie*ntrks) < 64) printf("iTrk = %i,  track (x,y,z,r)=(%.6f,%.6f,%.6f,%.6f) \n", it+ie*ntrks, x_,y_,z_,r_);
     }
   }

   avgpt = avgpt/float(nevts*ntrks);
   avgphi = avgphi/float(nevts*ntrks);
   avgtheta = avgtheta/float(nevts*ntrks);
   avgx = avgx/float(nevts*ntrks);
   avgy = avgy/float(nevts*ntrks);
   avgz = avgz/float(nevts*ntrks);
   avgr = avgr/float(nevts*ntrks);
   avgdx = avgdx/float(nevts*ntrks);
   avgdy = avgdy/float(nevts*ntrks);
   avgdz = avgdz/float(nevts*ntrks);
   avgdr = avgdr/float(nevts*ntrks);

   float stdx = 0, stdy = 0, stdz = 0, stdr = 0;
   float stddx = 0, stddy = 0, stddz = 0, stddr = 0;
   for (int ie=0;ie<nevts;++ie) {
     for (int it=0;it<ntrks;++it) {
       float x_ = x(outtrk,ie,it);
       float y_ = y(outtrk,ie,it);
       float z_ = z(outtrk,ie,it);
       float r_ = sqrtf(x_*x_ + y_*y_);
       float hx_ = inputhits[nlayer-1].pos[0];
       float hy_ = inputhits[nlayer-1].pos[1];
       float hz_ = inputhits[nlayer-1].pos[2];
       float hr_ = sqrtf(hx_*hx_ + hy_*hy_);
       if (std::isfinite(x_)==false ||
          std::isfinite(y_)==false ||
          std::isfinite(z_)==false
          ) {
        continue;
       }
       if (fabs( (x_-hx_)/hx_ )>1. ||
	   fabs( (y_-hy_)/hy_ )>1. ||
	   fabs( (z_-hz_)/hz_ )>1.) {
	 continue;
       }
       stdx += (x_-avgx)*(x_-avgx);
       stdy += (y_-avgy)*(y_-avgy);
       stdz += (z_-avgz)*(z_-avgz);
       stdr += (r_-avgr)*(r_-avgr);
       stddx += ((x_-hx_)/x_-avgdx)*((x_-hx_)/x_-avgdx);
       stddy += ((y_-hy_)/y_-avgdy)*((y_-hy_)/y_-avgdy);
       stddz += ((z_-hz_)/z_-avgdz)*((z_-hz_)/z_-avgdz);
       stddr += ((r_-hr_)/r_-avgdr)*((r_-hr_)/r_-avgdr);
     }
   }

   stdx = sqrtf(stdx/float(nevts*ntrks));
   stdy = sqrtf(stdy/float(nevts*ntrks));
   stdz = sqrtf(stdz/float(nevts*ntrks));
   stdr = sqrtf(stdr/float(nevts*ntrks));
   stddx = sqrtf(stddx/float(nevts*ntrks));
   stddy = sqrtf(stddy/float(nevts*ntrks));
   stddz = sqrtf(stddz/float(nevts*ntrks));
   stddr = sqrtf(stddr/float(nevts*ntrks));

   printf("track x avg=%.7f std/avg=%.7f\n", avgx, fabs(stdx/avgx));
   printf("track y avg=%.7f std/avg=%.7f\n", avgy, fabs(stdy/avgy));
   printf("track z avg=%.7f std/avg=%.7f\n", avgz, fabs(stdz/avgz));
   printf("track r avg=%.7f std/avg=%.7f\n", avgr, fabs(stdr/avgz));
   printf("track dx/x avg=%.7f std=%.7f\n", avgdx, stddx);
   printf("track dy/y avg=%.7f std=%.7f\n", avgdy, stddy);
   printf("track dz/z avg=%.7f std=%.7f\n", avgdz, stddz);
   printf("track dr/r avg=%.7f std=%.7f\n", avgdr, stddr);
   printf("track pt avg=%.7f\n", avgpt);
   printf("track phi avg=%.7f\n", avgphi);
   printf("track theta avg=%.7f\n", avgtheta);
   printf("number of tracks with nans=%i\n", nnans);
   printf("number of tracks failed=%i\n", nfail);

   return 0;
}
