#pragma once
#include <p2r_mpnx.h>

inline void MultHelixProp(const MP6x6F_ &a, const MP6x6SF_ &b, MP6x6F_ &c) {//ok

  c[ 0] = a[ 0]*b[ 0] + a[ 1]*b[ 1] + a[ 3]*b[ 6] + a[ 4]*b[10];
  c[ 1] = a[ 0]*b[ 1] + a[ 1]*b[ 2] + a[ 3]*b[ 7] + a[ 4]*b[11];
  c[ 2] = a[ 0]*b[ 3] + a[ 1]*b[ 4] + a[ 3]*b[ 8] + a[ 4]*b[12];
  c[ 3] = a[ 0]*b[ 6] + a[ 1]*b[ 7] + a[ 3]*b[ 9] + a[ 4]*b[13];
  c[ 4] = a[ 0]*b[10] + a[ 1]*b[11] + a[ 3]*b[13] + a[ 4]*b[14];
  c[ 5] = a[ 0]*b[15] + a[ 1]*b[16] + a[ 3]*b[18] + a[ 4]*b[19];
  c[ 6] = a[ 6]*b[ 0] + a[ 7]*b[ 1] + a[ 9]*b[ 6] + a[10]*b[10];
  c[ 7] = a[ 6]*b[ 1] + a[ 7]*b[ 2] + a[ 9]*b[ 7] + a[10]*b[11];
  c[ 8] = a[ 6]*b[ 3] + a[ 7]*b[ 4] + a[ 9]*b[ 8] + a[10]*b[12];
  c[ 9] = a[ 6]*b[ 6] + a[ 7]*b[ 7] + a[ 9]*b[ 9] + a[10]*b[13];
  c[10] = a[ 6]*b[10] + a[ 7]*b[11] + a[ 9]*b[13] + a[10]*b[14];
  c[11] = a[ 6]*b[15] + a[ 7]*b[16] + a[ 9]*b[18] + a[10]*b[19];
    
  c[12] = a[12]*b[ 0] + a[13]*b[ 1] + b[ 3] + a[15]*b[ 6] + a[16]*b[10] + a[17]*b[15];
  c[13] = a[12]*b[ 1] + a[13]*b[ 2] + b[ 4] + a[15]*b[ 7] + a[16]*b[11] + a[17]*b[16];
  c[14] = a[12]*b[ 3] + a[13]*b[ 4] + b[ 5] + a[15]*b[ 8] + a[16]*b[12] + a[17]*b[17];
  c[15] = a[12]*b[ 6] + a[13]*b[ 7] + b[ 8] + a[15]*b[ 9] + a[16]*b[13] + a[17]*b[18];
  c[16] = a[12]*b[10] + a[13]*b[11] + b[12] + a[15]*b[13] + a[16]*b[14] + a[17]*b[19];
  c[17] = a[12]*b[15] + a[13]*b[16] + b[17] + a[15]*b[18] + a[16]*b[19] + a[17]*b[20];
    
  c[18] = a[18]*b[ 0] + a[19]*b[ 1] + a[21]*b[ 6] + a[22]*b[10];
  c[19] = a[18]*b[ 1] + a[19]*b[ 2] + a[21]*b[ 7] + a[22]*b[11];
  c[20] = a[18]*b[ 3] + a[19]*b[ 4] + a[21]*b[ 8] + a[22]*b[12];
  c[21] = a[18]*b[ 6] + a[19]*b[ 7] + a[21]*b[ 9] + a[22]*b[13];
  c[22] = a[18]*b[10] + a[19]*b[11] + a[21]*b[13] + a[22]*b[14];
  c[23] = a[18]*b[15] + a[19]*b[16] + a[21]*b[18] + a[22]*b[19];
  c[24] = a[24]*b[ 0] + a[25]*b[ 1] + a[27]*b[ 6] + a[28]*b[10];
  c[25] = a[24]*b[ 1] + a[25]*b[ 2] + a[27]*b[ 7] + a[28]*b[11];
  c[26] = a[24]*b[ 3] + a[25]*b[ 4] + a[27]*b[ 8] + a[28]*b[12];
  c[27] = a[24]*b[ 6] + a[25]*b[ 7] + a[27]*b[ 9] + a[28]*b[13];
  c[28] = a[24]*b[10] + a[25]*b[11] + a[27]*b[13] + a[28]*b[14];
  c[29] = a[24]*b[15] + a[25]*b[16] + a[27]*b[18] + a[28]*b[19];
  c[30] = b[15];
  c[31] = b[16];
  c[32] = b[17];
  c[33] = b[18];
  c[34] = b[19];
  c[35] = b[20];    
  
  return;
}

inline void MultHelixPropTransp(const MP6x6F_ &a, const MP6x6F_ &b, MP6x6SF_ &c) {//

  c[ 0] = b[ 0]*a[ 0] + b[ 1]*a[ 1] + b[ 3]*a[ 3] + b[ 4]*a[ 4];
  c[ 1] = b[ 6]*a[ 0] + b[ 7]*a[ 1] + b[ 9]*a[ 3] + b[10]*a[ 4];
  c[ 2] = b[ 6]*a[ 6] + b[ 7]*a[ 7] + b[ 9]*a[ 9] + b[10]*a[10];
  c[ 3] = b[12]*a[ 0] + b[13]*a[ 1] + b[15]*a[ 3] + b[16]*a[ 4];
  c[ 4] = b[12]*a[ 6] + b[13]*a[ 7] + b[15]*a[ 9] + b[16]*a[10];
  c[ 5] = b[12]*a[12] + b[13]*a[13] + b[14] + b[15]*a[15] + b[16]*a[16] + b[17]*a[17];
  c[ 6] = b[18]*a[ 0] + b[19]*a[ 1] + b[21]*a[ 3] + b[22]*a[ 4];
  c[ 7] = b[18]*a[ 6] + b[19]*a[ 7] + b[21]*a[ 9] + b[22]*a[10];
  c[ 8] = b[18]*a[12] + b[19]*a[13] + b[20] + b[21]*a[15] + b[22]*a[16] + b[23]*a[17];
  c[ 9] = b[18]*a[18] + b[19]*a[19] + b[21]*a[21] + b[22]*a[22];
  c[10] = b[24]*a[ 0] + b[25]*a[ 1] + b[27]*a[ 3] + b[28]*a[ 4];
  c[11] = b[24]*a[ 6] + b[25]*a[ 7] + b[27]*a[ 9] + b[28]*a[10];
  c[12] = b[24]*a[12] + b[25]*a[13] + b[26] + b[27]*a[15] + b[28]*a[16] + b[29]*a[17];
  c[13] = b[24]*a[18] + b[25]*a[19] + b[27]*a[21] + b[28]*a[22];
  c[14] = b[24]*a[24] + b[25]*a[25] + b[27]*a[27] + b[28]*a[28];
  c[15] = b[30]*a[ 0] + b[31]*a[ 1] + b[33]*a[ 3] + b[34]*a[ 4];
  c[16] = b[30]*a[ 6] + b[31]*a[ 7] + b[33]*a[ 9] + b[34]*a[10];
  c[17] = b[30]*a[12] + b[31]*a[13] + b[32] + b[33]*a[15] + b[34]*a[16] + b[35]*a[17];
  c[18] = b[30]*a[18] + b[31]*a[19] + b[33]*a[21] + b[34]*a[22];
  c[19] = b[30]*a[24] + b[31]*a[25] + b[33]*a[27] + b[34]*a[28];
  c[20] = b[35];
  
  return;  
}

inline float hipo(const float x, const float y) {return std::sqrt(x*x + y*y);}

inline void KalmanUpdate(MP6x6SF_ &trkErr_, MP6F_ &inPar_, const MP3x3SF_ &hitErr_, const MP3F_ &msP_){	  
  
  MP1F_    rotT00;
  MP1F_    rotT01;
  MP2x2SF_ resErr_loc;
  //MP3x3SF_ resErr_glo;
  {   
    const auto msPX = msP_[iparX];
    const auto msPY = msP_[iparY];
    const auto inParX = inPar_[iparX];
    const auto inParY = inPar_[iparY];          
  
    const auto r = hipo(msPX, msPY);
    rotT00[0] = -(msPY + inParY) / (2*r);
    rotT01[0] =  (msPX + inParX) / (2*r);    
    
    resErr_loc[ 0] = (rotT00[0]*(trkErr_[0] + hitErr_[0]) +
                                    rotT01[0]*(trkErr_[1] + hitErr_[1]))*rotT00[0] +
                                   (rotT00[0]*(trkErr_[1] + hitErr_[1]) +
                                    rotT01[0]*(trkErr_[2] + hitErr_[2]))*rotT01[0];
    resErr_loc[ 1] = (trkErr_[3] + hitErr_[3])*rotT00[0] +
                                   (trkErr_[4] + hitErr_[4])*rotT01[0];
    resErr_loc[ 2] = (trkErr_[5] + hitErr_[5]);
  } 
  
  {
  
    const double det = (double)resErr_loc[0] * resErr_loc[2] -
                       (double)resErr_loc[1] * resErr_loc[1];
    const float s   = 1.f / det;
    const float tmp = s * resErr_loc[2];
    resErr_loc[1] *= -s;
    resErr_loc[2]  = s * resErr_loc[0];
    resErr_loc[0]  = tmp;  
  }     
  
  MP3x6_ kGain;
  
  {
    kGain[ 0] = trkErr_[ 0]*(rotT00[0]*resErr_loc[ 0]) +
	                        trkErr_[ 1]*(rotT01[0]*resErr_loc[ 0]) +
	                        trkErr_[ 3]*resErr_loc[ 1];
    kGain[ 1] = trkErr_[ 0]*(rotT00[0]*resErr_loc[ 1]) +
	                        trkErr_[ 1]*(rotT01[0]*resErr_loc[ 1]) +
	                        trkErr_[ 3]*resErr_loc[ 2];
    kGain[ 2] = 0;
    kGain[ 3] = trkErr_[ 1]*(rotT00[0]*resErr_loc[ 0]) +
	                        trkErr_[ 2]*(rotT01[0]*resErr_loc[ 0]) +
	                        trkErr_[ 4]*resErr_loc[ 1];
    kGain[ 4] = trkErr_[ 1]*(rotT00[0]*resErr_loc[ 1]) +
	                        trkErr_[ 2]*(rotT01[0]*resErr_loc[ 1]) +
	                        trkErr_[ 4]*resErr_loc[ 2];
    kGain[ 5] = 0;
    kGain[ 6] = trkErr_[ 3]*(rotT00[0]*resErr_loc[ 0]) +
	                        trkErr_[ 4]*(rotT01[0]*resErr_loc[ 0]) +
	                        trkErr_[ 5]*resErr_loc[ 1];
    kGain[ 7] = trkErr_[ 3]*(rotT00[0]*resErr_loc[ 1]) +
	                        trkErr_[ 4]*(rotT01[0]*resErr_loc[ 1]) +
	                        trkErr_[ 5]*resErr_loc[ 2];
    kGain[ 8] = 0;
    kGain[ 9] = trkErr_[ 6]*(rotT00[0]*resErr_loc[ 0]) +
	                        trkErr_[ 7]*(rotT01[0]*resErr_loc[ 0]) +
	                        trkErr_[ 8]*resErr_loc[ 1];
    kGain[10] = trkErr_[ 6]*(rotT00[0]*resErr_loc[ 1]) +
	                        trkErr_[ 7]*(rotT01[0]*resErr_loc[ 1]) +
	                        trkErr_[ 8]*resErr_loc[ 2];
    kGain[11] = 0;
    kGain[12] = trkErr_[10]*(rotT00[0]*resErr_loc[ 0]) +
	                        trkErr_[11]*(rotT01[0]*resErr_loc[ 0]) +
	                        trkErr_[12]*resErr_loc[ 1];
    kGain[13] = trkErr_[10]*(rotT00[0]*resErr_loc[ 1]) +
	                        trkErr_[11]*(rotT01[0]*resErr_loc[ 1]) +
	                        trkErr_[12]*resErr_loc[ 2];
    kGain[14] = 0;
    kGain[15] = trkErr_[15]*(rotT00[0]*resErr_loc[ 0]) +
	                        trkErr_[16]*(rotT01[0]*resErr_loc[ 0]) +
	                        trkErr_[17]*resErr_loc[ 1];
    kGain[16] = trkErr_[15]*(rotT00[0]*resErr_loc[ 1]) +
	                        trkErr_[16]*(rotT01[0]*resErr_loc[ 1]) +
	                        trkErr_[17]*resErr_loc[ 2];
    kGain[17] = 0;  
  }  
     
  MP2F_ res_loc;   
  {
    const auto msPX = msP_[iparX];
    const auto msPY = msP_[iparY];
    const auto msPZ = msP_[iparZ];    
    const auto inParX = inPar_[iparX];
    const auto inParY = inPar_[iparY];     
    const auto inParZ = inPar_[iparZ]; 
    
    const auto inParIpt   = inPar_[iparIpt];
    const auto inParPhi   = inPar_[iparPhi];
    const auto inParTheta = inPar_[iparTheta];            
    
    res_loc[0] =  rotT00[0]*(msPX - inParX) + rotT01[0]*(msPY - inParY);
    res_loc[1] =  msPZ - inParZ;

    inPar_[iparX]     = inParX + kGain[ 0] * res_loc[ 0] + kGain[ 1] * res_loc[ 1];
    inPar_[iparY]     = inParY + kGain[ 3] * res_loc[ 0] + kGain[ 4] * res_loc[ 1];
    inPar_[iparZ]     = inParZ + kGain[ 6] * res_loc[ 0] + kGain[ 7] * res_loc[ 1];
    inPar_[iparIpt]   = inParIpt + kGain[ 9] * res_loc[ 0] + kGain[10] * res_loc[ 1];
    inPar_[iparPhi]   = inParPhi + kGain[12] * res_loc[ 0] + kGain[13] * res_loc[ 1];
    inPar_[iparTheta] = inParTheta + kGain[15] * res_loc[ 0] + kGain[16] * res_loc[ 1];     
  }

  MP6x6SF_ newErr;
  {

     newErr[ 0] = kGain[ 0]*rotT00[0]*trkErr_[ 0] +
                         kGain[ 0]*rotT01[0]*trkErr_[ 1] +
                         kGain[ 1]*trkErr_[ 3];
     newErr[ 1] = kGain[ 3]*rotT00[0]*trkErr_[ 0] +
                         kGain[ 3]*rotT01[0]*trkErr_[ 1] +
                         kGain[ 4]*trkErr_[ 3];
     newErr[ 2] = kGain[ 3]*rotT00[0]*trkErr_[ 1] +
                         kGain[ 3]*rotT01[0]*trkErr_[ 2] +
                         kGain[ 4]*trkErr_[ 4];
     newErr[ 3] = kGain[ 6]*rotT00[0]*trkErr_[ 0] +
                         kGain[ 6]*rotT01[0]*trkErr_[ 1] +
                         kGain[ 7]*trkErr_[ 3];
     newErr[ 4] = kGain[ 6]*rotT00[0]*trkErr_[ 1] +
                         kGain[ 6]*rotT01[0]*trkErr_[ 2] +
                         kGain[ 7]*trkErr_[ 4];
     newErr[ 5] = kGain[ 6]*rotT00[0]*trkErr_[ 3] +
                         kGain[ 6]*rotT01[0]*trkErr_[ 4] +
                         kGain[ 7]*trkErr_[ 5];
     newErr[ 6] = kGain[ 9]*rotT00[0]*trkErr_[ 0] +
                         kGain[ 9]*rotT01[0]*trkErr_[ 1] +
                         kGain[10]*trkErr_[ 3];
     newErr[ 7] = kGain[ 9]*rotT00[0]*trkErr_[ 1] +
                         kGain[ 9]*rotT01[0]*trkErr_[ 2] +
                         kGain[10]*trkErr_[ 4];
     newErr[ 8] = kGain[ 9]*rotT00[0]*trkErr_[ 3] +
                         kGain[ 9]*rotT01[0]*trkErr_[ 4] +
                         kGain[10]*trkErr_[ 5];
     newErr[ 9] = kGain[ 9]*rotT00[0]*trkErr_[ 6] +
                         kGain[ 9]*rotT01[0]*trkErr_[ 7] +
                         kGain[10]*trkErr_[ 8];
     newErr[10] = kGain[12]*rotT00[0]*trkErr_[ 0] +
                         kGain[12]*rotT01[0]*trkErr_[ 1] +
                         kGain[13]*trkErr_[ 3];
     newErr[11] = kGain[12]*rotT00[0]*trkErr_[ 1] +
                         kGain[12]*rotT01[0]*trkErr_[ 2] +
                         kGain[13]*trkErr_[ 4];
     newErr[12] = kGain[12]*rotT00[0]*trkErr_[ 3] +
                         kGain[12]*rotT01[0]*trkErr_[ 4] +
                         kGain[13]*trkErr_[ 5];
     newErr[13] = kGain[12]*rotT00[0]*trkErr_[ 6] +
                         kGain[12]*rotT01[0]*trkErr_[ 7] +
                         kGain[13]*trkErr_[ 8];
     newErr[14] = kGain[12]*rotT00[0]*trkErr_[10] +
                         kGain[12]*rotT01[0]*trkErr_[11] +
                         kGain[13]*trkErr_[12];
     newErr[15] = kGain[15]*rotT00[0]*trkErr_[ 0] +
                         kGain[15]*rotT01[0]*trkErr_[ 1] +
                         kGain[16]*trkErr_[ 3];
     newErr[16] = kGain[15]*rotT00[0]*trkErr_[ 1] +
                         kGain[15]*rotT01[0]*trkErr_[ 2] +
                         kGain[16]*trkErr_[ 4];
     newErr[17] = kGain[15]*rotT00[0]*trkErr_[ 3] +
                         kGain[15]*rotT01[0]*trkErr_[ 4] +
                         kGain[16]*trkErr_[ 5];
     newErr[18] = kGain[15]*rotT00[0]*trkErr_[ 6] +
                         kGain[15]*rotT01[0]*trkErr_[ 7] +
                         kGain[16]*trkErr_[ 8];
     newErr[19] = kGain[15]*rotT00[0]*trkErr_[10] +
                         kGain[15]*rotT01[0]*trkErr_[11] +
                         kGain[16]*trkErr_[12];
     newErr[20] = kGain[15]*rotT00[0]*trkErr_[15] +
                         kGain[15]*rotT01[0]*trkErr_[16] +
                         kGain[16]*trkErr_[17];     
 #pragma unroll
     for (int i = 0; i < 21; i++){
       trkErr_[ i] = trkErr_[ i] - newErr[ i];
     }
   }
   //
   return;                 
}
                  

constexpr float kfact= 100/(-0.299792458*3.8112);
constexpr int Niter=5;

inline void propagateToR(const MP6x6SF_ &inErr_, const MP6F_ &inPar_, const MP1I_ &inChg_, 
                  const MP3F_ &msP_, MP6x6SF_ &outErr_, MP6F_ &outPar_) {
  //aux objects  
  MP6x6F_ errorProp;
  MP6x6F_ temp;
  
  auto PosInMtrx = [=] (int i, int j, int D) constexpr {return (i*D+j);};
  
  auto sincos4 = [] (const float x, float& sin, float& cos) {
    const float x2 = x*x;
    cos  = 1.f - 0.5f*x2 + 0.04166667f*x2*x2;
    sin  = x - 0.16666667f*x*x2;
  };
  
  {
    //initialize erroProp to identity matrix
    //for (int i=0;i<6;++i) errorProp.data[bsize*PosInMtrx(i,i,6) + it] = 1.f; 
    errorProp[PosInMtrx(0,0,6)] = 1.0f;
    errorProp[PosInMtrx(1,1,6)] = 1.0f;
    errorProp[PosInMtrx(2,2,6)] = 1.0f;
    errorProp[PosInMtrx(3,3,6)] = 1.0f;
    errorProp[PosInMtrx(4,4,6)] = 1.0f;
    errorProp[PosInMtrx(5,5,6)] = 1.0f;
    //
    const auto xin = inPar_[iparX];
    const auto yin = inPar_[iparY];     
    const auto zin = inPar_[iparZ]; 
    
    const auto iptin   = inPar_[iparIpt];
    const auto phiin   = inPar_[iparPhi];
    const auto thetain = inPar_[iparTheta]; 
    //
    auto r0 = hipo(xin, yin);
    const auto k = inChg_[0]*kfact;//?
    
    const auto xmsP = msP_[iparX];//?
    const auto ymsP = msP_[iparY];//?
    
    const auto r = hipo(xmsP, ymsP);    
    
    outPar_[iparX] = xin;
    outPar_[iparY] = yin;
    outPar_[iparZ] = zin;

    outPar_[iparIpt]   = iptin;
    outPar_[iparPhi]   = phiin;
    outPar_[iparTheta] = thetain;
 
    const auto kinv  = 1.f/k;
    const auto pt = 1.f/iptin;

    auto D = 0.f, cosa = 0.f, sina = 0.f, id = 0.f;
    //no trig approx here, phi can be large
    auto cosPorT = std::cos(phiin), sinPorT = std::sin(phiin);
    auto pxin = cosPorT*pt;
    auto pyin = sinPorT*pt;

    //derivatives initialized to value for first iteration, i.e. distance = r-r0in
    auto dDdx = r0 > 0.f ? -xin/r0 : 0.f;
    auto dDdy = r0 > 0.f ? -yin/r0 : 0.f;
    auto dDdipt = 0.;
    auto dDdphi = 0.;  
#pragma unroll    
    for (int i = 0; i < Niter; ++i)
    {
     //compute distance and path for the current iteration
      const auto xout = outPar_[iparX];
      const auto yout = outPar_[iparY];     
      
      r0 = hipo(xout, yout);
      id = (r-r0);
      D+=id;
      sincos4(id*iptin*kinv, sina, cosa);

      //update derivatives on total distance
      if (i+1 != Niter) {

	const auto oor0 = (r0>0.f && std::abs(r-r0)<0.0001f) ? 1.f/r0 : 0.f;

	const auto dadipt = id*kinv;

	const auto dadx = -xout*iptin*kinv*oor0;
	const auto dady = -yout*iptin*kinv*oor0;

	const auto pxca = pxin*cosa;
	const auto pxsa = pxin*sina;
	const auto pyca = pyin*cosa;
	const auto pysa = pyin*sina;

	auto tmp = k*dadx;
	dDdx   -= ( xout*(1.f + tmp*(pxca - pysa)) + yout*tmp*(pyca + pxsa) )*oor0;
	tmp = k*dady;
	dDdy   -= ( xout*tmp*(pxca - pysa) + yout*(1.f + tmp*(pyca + pxsa)) )*oor0;
	//now r0 depends on ipt and phi as well
	tmp = dadipt*iptin;
	dDdipt -= k*( xout*(pxca*tmp - pysa*tmp - pyca - pxsa + pyin) +
		      yout*(pyca*tmp + pxsa*tmp - pysa + pxca - pxin))*pt*oor0;
	dDdphi += k*( xout*(pysa - pxin + pxca) - yout*(pxsa - pyin + pyca))*oor0;
      } 
      
      //update parameters
      outPar_[iparX] = xout + k*(pxin*sina - pyin*(1.f-cosa));
      outPar_[iparY] = yout + k*(pyin*sina + pxin*(1.f-cosa));
      const float pxinold = pxin;//copy before overwriting
      pxin = pxin*cosa - pyin*sina;
      pyin = pyin*cosa + pxinold*sina;
  
    }
    //
    const auto alpha  = D*iptin*kinv;
    const auto dadx   = dDdx*iptin*kinv;
    const auto dady   = dDdy*iptin*kinv;
    const auto dadipt = (iptin*dDdipt + D)*kinv;
    const auto dadphi = dDdphi*iptin*kinv;

    sincos4(alpha, sina, cosa);
 
    errorProp[PosInMtrx(0,0,6)] = 1.f+k*dadx*(cosPorT*cosa-sinPorT*sina)*pt;
    errorProp[PosInMtrx(0,1,6)] =     k*dady*(cosPorT*cosa-sinPorT*sina)*pt;
    errorProp[PosInMtrx(0,2,6)] = 0.f;
    errorProp[PosInMtrx(0,3,6)] = k*(cosPorT*(iptin*dadipt*cosa-sina)+sinPorT*((1.f-cosa)-iptin*dadipt*sina))*pt*pt;
    errorProp[PosInMtrx(0,4,6)] = k*(cosPorT*dadphi*cosa - sinPorT*dadphi*sina - sinPorT*sina + cosPorT*cosa - cosPorT)*pt;
    errorProp[PosInMtrx(0,5,6)] = 0.f;

    errorProp[PosInMtrx(1,0,6)] =     k*dadx*(sinPorT*cosa+cosPorT*sina)*pt;
    errorProp[PosInMtrx(1,1,6)] = 1.f+k*dady*(sinPorT*cosa+cosPorT*sina)*pt;
    errorProp[PosInMtrx(1,2,6)] = 0.f;
    errorProp[PosInMtrx(1,3,6)] = k*(sinPorT*(iptin*dadipt*cosa-sina)+cosPorT*(iptin*dadipt*sina-(1.f-cosa)))*pt*pt;
    errorProp[PosInMtrx(1,4,6)] = k*(sinPorT*dadphi*cosa + cosPorT*dadphi*sina + sinPorT*cosa + cosPorT*sina - sinPorT)*pt;
    errorProp[PosInMtrx(1,5,6)] = 0.f;

    //no trig approx here, theta can be large
    cosPorT=std::cos(thetain);
    sinPorT=std::sin(thetain);
    //redefine sinPorT as 1./sinPorT to reduce the number of temporaries
    sinPorT = 1.f/sinPorT;

    outPar_[iparZ] = zin + k*alpha*cosPorT*pt*sinPorT;    

    errorProp[PosInMtrx(2,0,6)] = k*cosPorT*dadx*pt*sinPorT;
    errorProp[PosInMtrx(2,1,6)] = k*cosPorT*dady*pt*sinPorT;
    errorProp[PosInMtrx(2,2,6)] = 1.f;
    errorProp[PosInMtrx(2,3,6)] = k*cosPorT*(iptin*dadipt-alpha)*pt*pt*sinPorT;
    errorProp[PosInMtrx(2,4,6)] = k*dadphi*cosPorT*pt*sinPorT;
    errorProp[PosInMtrx(2,5,6)] =-k*alpha*pt*sinPorT*sinPorT;   
    //
    outPar_[iparIpt] = iptin;
 
    errorProp[PosInMtrx(3,0,6)] = 0.f;
    errorProp[PosInMtrx(3,1,6)] = 0.f;
    errorProp[PosInMtrx(3,2,6)] = 0.f;
    errorProp[PosInMtrx(3,3,6)] = 1.f;
    errorProp[PosInMtrx(3,4,6)] = 0.f;
    errorProp[PosInMtrx(3,5,6)] = 0.f; 
    
    outPar_[iparPhi] = phiin+alpha;
   
    errorProp[PosInMtrx(4,0,6)] = dadx;
    errorProp[PosInMtrx(4,1,6)] = dady;
    errorProp[PosInMtrx(4,2,6)] = 0.f;
    errorProp[PosInMtrx(4,3,6)] = dadipt;
    errorProp[PosInMtrx(4,4,6)] = 1.f+dadphi;
    errorProp[PosInMtrx(4,5,6)] = 0.f; 
  
    outPar_[iparTheta] = thetain;        

    errorProp[PosInMtrx(5,0,6)] = 0.f;
    errorProp[PosInMtrx(5,1,6)] = 0.f;
    errorProp[PosInMtrx(5,2,6)] = 0.f;
    errorProp[PosInMtrx(5,3,6)] = 0.f;
    errorProp[PosInMtrx(5,4,6)] = 0.f;
    errorProp[PosInMtrx(5,5,6)] = 1.f; 
                                 
  }
  
  MultHelixProp(errorProp, inErr_, temp);
  MultHelixPropTransp(errorProp, temp, outErr_);  
  
  return;
}



