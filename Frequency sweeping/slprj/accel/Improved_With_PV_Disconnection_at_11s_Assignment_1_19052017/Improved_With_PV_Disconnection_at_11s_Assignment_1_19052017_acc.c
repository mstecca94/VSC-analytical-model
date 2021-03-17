

    
          /*
     * This file is not available for use in any application other than as a
     * MATLAB(R) MEX file for use with the Simulink(R) product.
     */

  
    /*
    * Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc.c
    * 
      * Academic License - for use in teaching, academic research, and meeting
* course requirements at degree granting institutions only.  Not for
* government, commercial, or other organizational use. 
  * 
  * Code generation for model "Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc".
  *
  * Model version              : 1.442
  * Simulink Coder version : 9.1 (R2019a) 23-Nov-2018
  * C source code generated on : Thu May 14 09:53:39 2020
 * 
 * Target selection: accel.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Emulation hardware selection: 
 *    Differs from embedded hardware (MATLAB Host)
 * Code generation objectives: Unspecified
 * Validation result: Not run
    */



    

        #include <math.h>

        #include "Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc.h"


        #include "Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_private.h"


#include <stdio.h>
#include "slexec_vm_simstruct_bridge.h"
#include "slexec_vm_zc_functions.h"
#include "slexec_vm_lookup_functions.h"
#include "slsv_diagnostic_codegen_c_api.h"
#include "simtarget/slSimTgtMdlrefSfcnBridge.h"
      #include "simstruc.h"
      #include "fixedpoint.h"


  

  

  

  #define CodeFormat S-Function
#define AccDefine1 Accelerator_S-Function


  

  

  

    #include "simtarget/slAccSfcnBridge.h"


  

  

  

  

  

  

  

    
                  
      #ifndef __RTW_UTFREE__  
      extern void * utMalloc(size_t);
      extern void   utFree(void *); 
      #endif
      
      
  
    /* Buffer management routine for variable delay block */
    boolean_T  Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf(
int_T       *bufSzPtr,        /* in/out - circular buffer size                 */
int_T       *tailPtr,         /* in/out - tail of circular buffer              */
int_T       *headPtr,         /* in/out - head of circular buffer              */
int_T       *lastPtr,         /* in/out - same logical 'last' referenced index */
real_T      tMinusDelay,      /* in     - last point we are looking at   */
real_T      **tBufPtr,        /* in/out - larger buffer for time         */
real_T      **uBufPtr,        /* in/out - larger buffer for input        */
real_T      **xBufPtr,        /* in/out - larger buffer for state        */
boolean_T   isfixedbuf,       /* in     - fixed buffer size enable       */
boolean_T istransportdelay,   /* in     - block acts as transport dela y */
int_T     *maxNewBufSzPtr)
    {
      
      int_T  testIdx;
      int_T  tail  = *tailPtr;
      int_T  bufSz = *bufSzPtr;
      real_T *tBuf = *tBufPtr;
      real_T *xBuf = (NULL); 
      
      int_T    numBuffer = 2;
      if (istransportdelay){
        numBuffer =3 ;
        xBuf= *xBufPtr;
      }
      
      /*    Get testIdx, the index of the second oldest data point and
      *    see if this is older than current sim time minus applied delay,
      *    used to see if we can move tail forward 
      */
      testIdx = (tail < (bufSz - 1)) ? (tail + 1) : 0;
      
      if ( (tMinusDelay <= tBuf[testIdx]) && !isfixedbuf) {
        int_T  j;
        real_T *tempT;
        real_T *tempU;
        real_T *tempX = (NULL);    
        
        real_T *uBuf     = *uBufPtr;
        int_T  newBufSz  = bufSz + 1024;
        
        if (newBufSz > *maxNewBufSzPtr) {
          *maxNewBufSzPtr = newBufSz; /* save for warning*/
        }
        
        tempU = (real_T*)utMalloc(numBuffer*newBufSz*sizeof(real_T));        
        
        if (tempU == (NULL)){
          return (false);
        }     
        tempT = tempU + newBufSz;
        if(istransportdelay) tempX = tempT + newBufSz;
        
        for (j = tail; j < bufSz; j++) {
          tempT[j - tail] = tBuf[j];
          tempU[j - tail] = uBuf[j];
          if (istransportdelay)  
          tempX[j - tail] = xBuf[j];
        }
        for (j = 0; j < tail; j++) {
          tempT[j + bufSz - tail] = tBuf[j];
          tempU[j + bufSz - tail] = uBuf[j];
          if (istransportdelay)  
          tempX[j + bufSz - tail] = xBuf[j];
        }
        
        if (*lastPtr> tail) 
        {
          *lastPtr -= tail;
        } else {
          *lastPtr += (bufSz - tail);
        }
        *tailPtr= 0;
        *headPtr = bufSz;
        
        utFree(uBuf);
        
        *bufSzPtr = newBufSz;
        *tBufPtr  = tempT;
        *uBufPtr  = tempU;
        if (istransportdelay) *xBufPtr  = tempX;
        
      }else {
        *tailPtr = testIdx; /* move tail forward */
      }
      
      return(true);      


    }
        
            /* 
 * Time delay interpolation routine
 * 
 * The linear interpolation is performed using the formula:
 * 
 *          (t2 - tMinusDelay)         (tMinusDelay - t1)
 * u(t)  =  ----------------- * u1  +  ------------------- * u2
 *              (t2 - t1)                  (t2 - t1)
 */

        real_T Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
      real_T     tMinusDelay,           /* tMinusDelay = currentSimTime - delay */
      real_T     tStart,
      real_T     *tBuf,
      real_T     *uBuf,
      int_T      bufSz,
      int_T      *lastIdx,
      int_T      oldestIdx,
      int_T      newIdx,
      real_T     initOutput,
      boolean_T  discrete,
      boolean_T  minorStepAndTAtLastMajorOutput)

    {
      
      int_T i;
      real_T yout, t1, t2, u1, u2;

      /*
      * If there is only one data point in the buffer, this data point must be
      * the t= 0 and tMinusDelay > t0, it ask for something unknown. The best
      * guess if initial output as well
      */
      
      if ( (newIdx == 0) && (oldestIdx ==0 ) && (tMinusDelay > tStart))  return initOutput;
      
      /*
      * If tMinusDelay is less than zero, should output initial value
      */
      
      if(tMinusDelay <= tStart) return initOutput;
      
      /* For fixed buffer extrapolation: 
      * if tMinusDelay is small than the time at oldestIdx, if discrete, output
      * tailptr value,  else use tailptr and tailptr+1 value to extrapolate
      * It is also for fixed buffer. Note: The same condition can happen for transport delay block where 
      * use tStart and and t[tail] other than using t[tail] and t[tail+1]. 
      * See below
      */
      
      if ( (tMinusDelay <= tBuf[oldestIdx] ) ) {
        if (discrete){
          return(uBuf[oldestIdx]);
        }else{
          int_T tempIdx=oldestIdx + 1;
          if(oldestIdx == bufSz-1) tempIdx = 0;
          t1=tBuf[oldestIdx];
          t2=tBuf[tempIdx];
          u1=uBuf[oldestIdx];
          u2=uBuf[tempIdx];
          
          if (t2 == t1) {
            if (tMinusDelay >= t2) {
              yout = u2;
            }else {
              yout = u1;
            }
          }
          else{
            real_T f1 = (t2-tMinusDelay) / (t2-t1);
            real_T f2 = 1.0 - f1;
            
            /*
            * Use Lagrange's interpolation formula.  Exact outputs at t1, t2.
            */
            yout = f1*u1 + f2*u2;
          }
          return yout;
        }
      }
        
      /*
      * When block does not have direct feedthrough, we use the table of
      * values to extrapolate off the end of the table for delays that are less
      * than 0 (less then step size).  This is not completely accurate.  The
      * chain of events is as follows for a given time t.  Major output - look 
      * in table.  Update - add entry to table.  Now, if we call the output at 
      * time t again, there is a new entry in the table. For very small delays,
      * this means that we will have a different answer from the previous call 
      * to the output fcn at the same time t.  The following code prevents this 
      * from happening.
      */
        
      if (minorStepAndTAtLastMajorOutput) {
        /* pretend that the new entry has not been added to table */
        if (newIdx != 0) {
          if (*lastIdx == newIdx) {
            (*lastIdx)--;
          }
          newIdx--;
        } else {
          if (*lastIdx == newIdx) {
            *lastIdx = bufSz-1;
          }
          newIdx = bufSz - 1;
        }
      }
          
      i = *lastIdx;
      if (tBuf[i] < tMinusDelay) {
        /* Look forward starting at last index */
        while (tBuf[i] < tMinusDelay) {
          
          /* May occur if the delay is less than step-size - extrapolate */
          if (i == newIdx) break;
          
          i = ( i < (bufSz-1) ) ? (i+1) : 0; /* move through buffer */
              
        }
      } else {
        /*
        * Look backwards starting at last index which can happen when the
        * delay time increases.
        */
        while (tBuf[i] >= tMinusDelay) {
          
          /*
          * Due to the entry condition at top of function, we
          * should never hit the end.
          */
          
          i = (i > 0) ? i-1 : (bufSz-1); /* move through buffer */
          
        }
        i = ( i < (bufSz-1) ) ? (i+1) : 0;
      }
          
      *lastIdx = i;
      
      if (discrete) {
        /*
        * tempEps = 128 * eps;
        * localEps = max(tempEps, tempEps*fabs(tBuf[i]))/2;
        */
        double tempEps  = (DBL_EPSILON) * 128.0;
        double localEps = tempEps * muDoubleScalarAbs(tBuf[i]);
        if (tempEps > localEps) {
          localEps = tempEps;
        }
        localEps = localEps / 2.0;
        
        if (tMinusDelay >= (tBuf[i] - localEps)) {
          yout = uBuf[i];
        } else {
          if (i == 0) {
            yout = uBuf[bufSz-1];
          } else {
            yout = uBuf[i-1];
          }
        }
      } else {
        if (i == 0) {
          t1 = tBuf[bufSz-1];
          u1 = uBuf[bufSz-1];
        } else {
          t1 = tBuf[i-1];
          u1 = uBuf[i-1];
        }
        
        t2 = tBuf[i];
        u2 = uBuf[i];
        
        if (t2 == t1) {
          if (tMinusDelay >= t2) {
            yout = u2;
          } else {
            yout = u1;
          }
        } else {
          real_T f1 = (t2-tMinusDelay) / (t2-t1);
          real_T f2 = 1.0 - f1;
          
          /*
          * Use Lagrange's interpolation formula.  Exact outputs at t1, t2.
          */
          yout = f1*u1 + f2*u2;
        }
      }
         
      return(yout); 
          
    }
        
          void rt_ssGetBlockPath(SimStruct* S, int_T sysIdx, int_T blkIdx, char_T **path) { _ssGetBlockPath(S, sysIdx, blkIdx, path); }
  void rt_ssSet_slErrMsg(SimStruct* S, void* diag) { if(!_ssIsErrorStatusAslErrMsg(S)) {_ssSet_slErrMsg(S, diag);} else {_ssDiscardDiagnostic(S, diag);}}
  void rt_ssReportDiagnosticAsWarning(SimStruct* S, void* diag) { _ssReportDiagnosticAsWarning(S, diag); }

    

    

  
  
  
                                        
                
    
        
        
     

      
         void rt_invd5x5_snf(const real_T u[25], real_T y[25])
    {
          int8_T p[5];
real_T A[25];
int8_T ipiv[5];
int32_T jj;
int32_T jp1j;
int32_T idxmax;
int32_T ix;
real_T smax;
real_T s;
int32_T iy;
int32_T jA;
int32_T ijA;
int32_T pipk;
int32_T jBcol;
for (iy = 0; iy < 25; iy++) {
    y[iy] = 0.0;
    A[iy] = u[iy];
}
for (iy = 0; iy < 5; iy++) {
    ipiv[iy] = (int8_T)(1 + iy);
}
for (pipk = 0; pipk < 4; pipk++) {
    jBcol = pipk * 6 + 1;
    jj = pipk * 6;
    jp1j = jBcol + 1;
    iy = 5 - pipk;
    idxmax = 1;
    ix = jBcol - 1;
    smax = muDoubleScalarAbs(A[jj]);
    for (jA = 2; jA <= iy; jA++) {
        ix++;
        s = muDoubleScalarAbs(A[ix]);
        if (s > smax) {
            idxmax = jA;
            smax = s;
        }
    }
    if (A[(jBcol + idxmax) - 2] != 0.0) {
        if (idxmax - 1 != 0) {
            ipiv[pipk] = (int8_T)(pipk + idxmax);
            ix = pipk + 1;
            iy = pipk + idxmax;
            for (jA = 0; jA < 5; jA++) {
                smax = A[ix - 1];
                A[ix - 1] = A[iy - 1];
                A[iy - 1] = smax;
                ix += 5;
                iy += 5;
            }
        }
        iy = jp1j - pipk;
        for (ix = jp1j; ix <= iy + 3; ix++) {
            A[ix - 1] /= A[jj];
        }
    }
    jBcol = 3 - pipk;
    jA = jj;
    jj += 5;
    for (idxmax = 0; idxmax <= jBcol; idxmax++) {
        if (A[jj] != 0.0) {
            smax = -A[jj];
            ix = jp1j;
            iy = jA - pipk;
            for (ijA = jA + 7; ijA <= iy + 10; ijA++) {
                A[ijA - 1] += A[ix - 1] * smax;
                ix++;
            }
        }
        jj += 5;
        jA += 5;
    }
}
for (iy = 0; iy < 5; iy++) {
    p[iy] = (int8_T)(1 + iy);
}
if (ipiv[0] > 1) {
    pipk = p[ipiv[0] - 1];
    p[ipiv[0] - 1] = p[0];
    p[0] = (int8_T)pipk;
}
if (ipiv[1] > 2) {
    pipk = p[ipiv[1] - 1];
    p[ipiv[1] - 1] = p[1];
    p[1] = (int8_T)pipk;
}
if (ipiv[2] > 3) {
    pipk = p[ipiv[2] - 1];
    p[ipiv[2] - 1] = p[2];
    p[2] = (int8_T)pipk;
}
if (ipiv[3] > 4) {
    pipk = p[4];
    p[4] = p[3];
    p[3] = (int8_T)pipk;
}
for (jA = 0; jA < 5; jA++) {
    jj = 1 + jA;
    jBcol = p[jj - 1];
    y[(jj + 5 * (p[jj - 1] - 1)) - 1] = 1.0;
    for (pipk = jj; pipk < 6; pipk++) {
        if (y[((jBcol - 1) * 5 + pipk) - 1] != 0.0) {
            for (ix = pipk + 1; ix < 6; ix++) {
                y[(ix + 5 * (jBcol - 1)) - 1] -= y[((jBcol - 1) * 5 + pipk) - 1] * A[((pipk - 1) * 5 + ix) - 1];
            }
        }
    }
}
for (pipk = 0; pipk < 5; pipk++) {
    jBcol = 5 * pipk;
    for (jA = 4; jA >= 0; jA--) {
        jj = 5 * jA;
        if (y[jA + jBcol] != 0.0) {
            y[jA + jBcol] /= A[jA + jj];
            iy = jA - 1;
            for (ix = 0; ix <= iy; ix++) {
                jp1j = 1 + ix;
                y[(jp1j + jBcol) - 1] -= A[(jp1j + jj) - 1] * y[jA + jBcol];
            }
        }
    }
}


    }
      

      

      
      
	      /* 
 * System initialize for enable system:
 *    '<S73>/Saturation'
 *    '<S835>/Saturation'
 */

          
  
  
      
  
      void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_Saturation_Init(SimStruct *S, DW_Saturation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localDW, P_Saturation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localP)
  {
  
        
  
    
  

                
  

                    
                localDW->Lmd_sat_DSTATE = localP->P_9;
localDW->Lmq_sat_DSTATE = localP->P_2;






    
  

          }
      
  



       
    
  

          /* 
 * Outputs for enable system:
 *    '<S73>/Saturation'
 *    '<S835>/Saturation'
 */

      
            
  
  
      
  
      void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_Saturation(SimStruct *S, real_T rtu_Enable, const real_T rtu_phi[5], B_Saturation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localB, DW_Saturation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localDW, P_Saturation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localP)
  {
  
                
    
real_T rtb_B_3_24_0[25];
real_T B_3_24_0[25];
int32_T i;
int32_T i_0;
real_T B_3_4_0_idx_2;
real_T B_3_4_0_idx_0;
real_T B_3_4_0_idx_1;
int32_T isHit;

    
  

                    
  
    
  
    
  

    
                   
                  


if (rtu_Enable > 0.0) {
    isHit = ssIsSampleHit(S, 2, 0);
    if (isHit != 0) {
        localB->B_3_3_0 = localB->B_3_0_0;
        B_3_4_0_idx_0 = localP->P_7[0] * rtu_phi[1];
        B_3_4_0_idx_1 = localP->P_7[1] * rtu_phi[2];
        B_3_4_0_idx_2 = localP->P_7[2] * rtu_phi[3];
        localB->B_3_10_0 = 1.0 / ((localB->B_3_6_0[0] + localB->B_3_6_0[1]) + 1.0 / localDW->Lmd_sat_DSTATE);
        B_3_4_0_idx_0 = muDoubleScalarAbs(((B_3_4_0_idx_0 + B_3_4_0_idx_1) + B_3_4_0_idx_2) * localB->B_3_10_0);
        B_3_4_0_idx_1 = look1_pbinlxpw(B_3_4_0_idx_0, localP->P_11, localP->P_10, &localDW->m_bpIndex, 1U);
        if (B_3_4_0_idx_1 != 0.0) {
            B_3_4_0_idx_0 /= B_3_4_0_idx_1;
        } else {
            B_3_4_0_idx_0 = localB->B_3_14_0;
        }
        localB->B_3_17_0 = localP->P_13 * B_3_4_0_idx_0;
        if (localB->B_3_3_0) {
            isHit = ssIsSampleHit(S, 2, 0);
            if (isHit != 0) {
                B_3_4_0_idx_0 = localP->P_0[0] * rtu_phi[0];
                B_3_4_0_idx_1 = localP->P_0[1] * rtu_phi[4];
                localB->B_2_6_0 = 1.0 / (((localB->B_2_2_0[0] + localB->B_2_2_0[1]) + localB->B_2_2_0[2]) + 1.0 / localDW->Lmq_sat_DSTATE);
                B_3_4_0_idx_0 = muDoubleScalarAbs((B_3_4_0_idx_0 + B_3_4_0_idx_1) * localB->B_2_6_0);
                B_3_4_0_idx_1 = look1_pbinlxpw(B_3_4_0_idx_0, localP->P_4, localP->P_3, &localDW->m_bpIndex_c, 1U);
                if (B_3_4_0_idx_1 != 0.0) {
                    B_3_4_0_idx_0 /= B_3_4_0_idx_1;
                } else {
                    B_3_4_0_idx_0 = localB->B_2_10_0;
                }
                localB->B_2_13_0 = localP->P_6 * B_3_4_0_idx_0;
            }
            srUpdateBC(localDW->Lmq_sat_SubsysRanBC);
        }
        if (localB->B_3_1_0) {
            localB->B_3_22_0 = localB->B_2_13_0;
        } else {
            localB->B_3_22_0 = localB->B_3_21_0;
        }
        memcpy(&rtb_B_3_24_0[0], &localB->B_3_20_0[0], 25U * sizeof(real_T));
        rtb_B_3_24_0[0] = localB->B_3_22_0;
        rtb_B_3_24_0[4] = localB->B_3_22_0;
        rtb_B_3_24_0[20] = localB->B_3_22_0;
        rtb_B_3_24_0[24] = localB->B_3_22_0;
        for (isHit = 0; isHit < 3; isHit++) {
            rtb_B_3_24_0[1 + 5 * (1 + isHit)] = localB->B_3_17_0;
            rtb_B_3_24_0[2 + 5 * (1 + isHit)] = localB->B_3_17_0;
            rtb_B_3_24_0[3 + 5 * (1 + isHit)] = localB->B_3_17_0;
        }
        for (isHit = 0; isHit < 25; isHit++) {
            B_3_24_0[isHit] = rtb_B_3_24_0[isHit] + localB->B_3_25_0[isHit];
        }
        rt_invd5x5_snf(B_3_24_0, localB->B_3_27_0);
        for (isHit = 0; isHit < 5; isHit++) {
            for (i = 0; i < 5; i++) {
                localB->B_3_28_0[i + 5 * isHit] = 0.0;
                for (i_0 = 0; i_0 < 5; i_0++) {
                    localB->B_3_28_0[i + 5 * isHit] += localB->B_3_19_0[5 * i_0 + i] * localB->B_3_27_0[5 * isHit + i_0];
                }
            }
        }
        if (localB->B_3_2_0) {
            localB->B_3_30_0 = localB->B_2_6_0;
        } else {
            localB->B_3_30_0 = localB->B_3_29_0;
        }
    }
    if (ssIsMajorTimeStep(S) != 0) {
        srUpdateBC(localDW->Saturation_SubsysRanBC);
    }
}






    
  

        
  
  
    
        }
      
  






    
  

          /* 
 * Outputs for enable system:
 *    '<S73>/Saturation'
 *    '<S835>/Saturation'
 */

      
            
  
  
      
  
      void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_SaturationTID5(SimStruct *S, B_Saturation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localB, P_Saturation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localP)
  {
  
                
    
    
  

                
  
    
  
    
  


                   
                  localB->B_3_0_0 = localP->P_19;
localB->B_3_1_0 = localP->P_20;
localB->B_3_2_0 = localP->P_21;
localB->B_3_6_0[0] = localP->P_8[0];
localB->B_3_6_0[1] = localP->P_8[1];
localB->B_3_14_0 = localP->P_12;
localB->B_2_2_0[0] = localP->P_1[0];
localB->B_2_2_0[1] = localP->P_1[1];
localB->B_2_2_0[2] = localP->P_1[2];
localB->B_2_10_0 = localP->P_5;
localB->B_3_21_0 = localP->P_16;
memcpy(&localB->B_3_19_0[0], &localP->P_14[0], 25U * sizeof(real_T));
memcpy(&localB->B_3_20_0[0], &localP->P_15[0], 25U * sizeof(real_T));
memcpy(&localB->B_3_25_0[0], &localP->P_17[0], 25U * sizeof(real_T));
localB->B_3_29_0 = localP->P_18;






    
  

    
  
  
    
        }
      
  








       
    
  

          /* 
 * Update for enable system:
 *    '<S73>/Saturation'
 *    '<S835>/Saturation'
 */

      
            
  
  
      
  
      void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_Saturation_Update(SimStruct *S, real_T rtu_Enable, B_Saturation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localB, DW_Saturation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localDW)
  {
  
      
        
    
int32_T isHit;
  
    
  

                
      
    
  

                      
                  if (rtu_Enable > 0.0) {
    isHit = ssIsSampleHit(S, 2, 0);
    if (isHit != 0) {
        localDW->Lmd_sat_DSTATE = localB->B_3_17_0;
        if (localB->B_3_3_0) {
            isHit = ssIsSampleHit(S, 2, 0);
            if (isHit != 0) {
                localDW->Lmq_sat_DSTATE = localB->B_2_13_0;
            }
        }
    }
}






    
  

      
  

            }
      
  






    
  






            /* 
 * Termination for enable system:
 *    '<S73>/Saturation'
 *    '<S835>/Saturation'
 */

          
  
  
      
  
      void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_Saturation_Term(SimStruct *const S)
  {
  
          
    
    
  

                
  

                  
  

          }
      
  



      
                                                                                                                                                                                                                          
      
	      /* 
 * Enable for action system:
 *    '<S1>/If Action Subsystem2'
 *    '<S1>/If Action Subsystem3'
 */

          
  
  
      
  
      void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_IfActionSubsystem2_Enable(SimStruct *S)
  {
  
        
  
    
  

                
  

                    
                if (ssGetTaskTime(S, 0) != ssGetTStart(S)) {
    ssSetBlockStateForSolverChangedAtMajorStep(S);
}
  





    
  

          }
      
  



	      /* 
 * Disable for action system:
 *    '<S1>/If Action Subsystem2'
 *    '<S1>/If Action Subsystem3'
 */

          
  
  
      
  
      void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_IfActionSubsystem2_Disable(SimStruct *S)
  {
  
        
  
    
  

                    
  

                    
                ssSetBlockStateForSolverChangedAtMajorStep(S);
  





    
  

          }
      
  



       
    
  

          /* 
 * Output and update for action system:
 *    '<S1>/If Action Subsystem2'
 *    '<S1>/If Action Subsystem3'
 */

      
              
  
  
      
  
      void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_IfActionSubsystem2(SimStruct *S, real_T rtu_In1, real_T *rty_Out1)
  {
  
                      
      
int32_T isHit;

      
      
  

                          
      
  

                      
                isHit = ssIsSampleHit(S, 1, 0);
if (isHit != 0) {
    *rty_Out1 = rtu_In1;
}






      
  

        
    
      
          }
      
  



                





            /* 
 * Termination for action system:
 *    '<S1>/If Action Subsystem2'
 *    '<S1>/If Action Subsystem3'
 */

          
  
  
      
  
      void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_IfActionSubsystem2_Term(SimStruct *const S)
  {
  
          
    
    
  

                
  

                  
  

          }
      
  



      
                    
      
	      /* 
 * Disable for enable system:
 *    '<S659>/Neg. Seq. Computation'
 *    '<S659>/Pos. Seq. Computation'
 *    '<S701>/Neg. Seq. Computation'
 *    '<S701>/Pos. Seq. Computation'
 *    '<S743>/Neg. Seq. Computation'
 *    '<S743>/Pos. Seq. Computation'
 *    '<S785>/Neg. Seq. Computation'
 *    '<S785>/Pos. Seq. Computation'
 */

          
  
  
      
  
      void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_NegSeqComputation_Disable(SimStruct *S, DW_NegSeqComputation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localDW)
  {
  
        
  
    
  

                    
  

                    
                localDW->NegSeqComputation_MODE = false;





    
  

          }
      
  



       
    
  

          /* 
 * Output and update for enable system:
 *    '<S659>/Neg. Seq. Computation'
 *    '<S659>/Pos. Seq. Computation'
 *    '<S701>/Neg. Seq. Computation'
 *    '<S701>/Pos. Seq. Computation'
 *    '<S743>/Neg. Seq. Computation'
 *    '<S743>/Pos. Seq. Computation'
 *    '<S785>/Neg. Seq. Computation'
 *    '<S785>/Pos. Seq. Computation'
 */

      
              
  
  
      
  
      void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_NegSeqComputation(SimStruct *S, real_T rtu_Enable, creal_T rtu_In, creal_T rtu_In_h, creal_T rtu_In_l, B_NegSeqComputation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localB, DW_NegSeqComputation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localDW, P_NegSeqComputation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localP)
  {
  
                      
      
int32_T isHit;

      
      
  

                          
      
  

                      
                
isHit = ssIsSampleHit(S, 1, 0);
if ((isHit != 0) && (ssIsMajorTimeStep(S) != 0)) {
    if (rtu_Enable > 0.0) {
        if (!localDW->NegSeqComputation_MODE) {
            if (ssGetTaskTime(S, 1) != ssGetTStart(S)) {
                ssSetBlockStateForSolverChangedAtMajorStep(S);
            }
  
            localDW->NegSeqComputation_MODE = true;
        }
    } else {
        if (localDW->NegSeqComputation_MODE) {
            ssSetBlockStateForSolverChangedAtMajorStep(S);
  



      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_NegSeqComputation_Disable(S, localDW);      

        }
    }
}
if (localDW->NegSeqComputation_MODE) {
    localB->B_39_0_0[0].re = localP->P_1[0].re * rtu_In.re - localP->P_1[0].im * rtu_In.im;
    localB->B_39_0_0[0].im = localP->P_1[0].re * rtu_In.im + localP->P_1[0].im * rtu_In.re;
    localB->B_39_0_0[1].re = localP->P_1[1].re * rtu_In_h.re - localP->P_1[1].im * rtu_In_h.im;
    localB->B_39_0_0[1].im = localP->P_1[1].re * rtu_In_h.im + localP->P_1[1].im * rtu_In_h.re;
    localB->B_39_0_0[2].re = localP->P_1[2].re * rtu_In_l.re - localP->P_1[2].im * rtu_In_l.im;
    localB->B_39_0_0[2].im = localP->P_1[2].re * rtu_In_l.im + localP->P_1[2].im * rtu_In_l.re;
    localB->B_39_1_0.re = (localB->B_39_0_0[0].re + localB->B_39_0_0[1].re) + localB->B_39_0_0[2].re;
    localB->B_39_1_0.im = (localB->B_39_0_0[0].im + localB->B_39_0_0[1].im) + localB->B_39_0_0[2].im;
    localB->B_39_2_0.re = localP->P_0 * localB->B_39_1_0.re;
    localB->B_39_2_0.im = localP->P_0 * localB->B_39_1_0.im;
    if (ssIsMajorTimeStep(S) != 0) {
        srUpdateBC(localDW->NegSeqComputation_SubsysRanBC);
    }
}






      
  

        
    
      
          }
      
  



                





            /* 
 * Termination for enable system:
 *    '<S659>/Neg. Seq. Computation'
 *    '<S659>/Pos. Seq. Computation'
 *    '<S701>/Neg. Seq. Computation'
 *    '<S701>/Pos. Seq. Computation'
 *    '<S743>/Neg. Seq. Computation'
 *    '<S743>/Pos. Seq. Computation'
 *    '<S785>/Neg. Seq. Computation'
 *    '<S785>/Pos. Seq. Computation'
 */

          
  
  
      
  
      void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_NegSeqComputation_Term(SimStruct *const S)
  {
  
          
    
    
  

                
  

                  
  

          }
      
  



      
                    
      
	      /* 
 * Disable for enable system:
 *    '<S659>/Zero Seq. Computation'
 *    '<S701>/Zero Seq. Computation'
 *    '<S743>/Zero Seq. Computation'
 *    '<S785>/Zero Seq. Computation'
 */

          
  
  
      
  
      void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_ZeroSeqComputation_Disable(SimStruct *S, DW_ZeroSeqComputation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localDW)
  {
  
        
  
    
  

                    
  

                    
                localDW->ZeroSeqComputation_MODE = false;





    
  

          }
      
  



       
    
  

          /* 
 * Output and update for enable system:
 *    '<S659>/Zero Seq. Computation'
 *    '<S701>/Zero Seq. Computation'
 *    '<S743>/Zero Seq. Computation'
 *    '<S785>/Zero Seq. Computation'
 */

      
              
  
  
      
  
      void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_ZeroSeqComputation(SimStruct *S, real_T rtu_Enable, creal_T rtu_In, creal_T rtu_In_f, creal_T rtu_In_o, B_ZeroSeqComputation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localB, DW_ZeroSeqComputation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localDW, P_ZeroSeqComputation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localP)
  {
  
                      
      
int32_T isHit;

      
      
  

                          
      
  

                      
                
isHit = ssIsSampleHit(S, 1, 0);
if ((isHit != 0) && (ssIsMajorTimeStep(S) != 0)) {
    if (rtu_Enable > 0.0) {
        if (!localDW->ZeroSeqComputation_MODE) {
            if (ssGetTaskTime(S, 1) != ssGetTStart(S)) {
                ssSetBlockStateForSolverChangedAtMajorStep(S);
            }
  
            localDW->ZeroSeqComputation_MODE = true;
        }
    } else {
        if (localDW->ZeroSeqComputation_MODE) {
            ssSetBlockStateForSolverChangedAtMajorStep(S);
  



      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_ZeroSeqComputation_Disable(S, localDW);      

        }
    }
}
if (localDW->ZeroSeqComputation_MODE) {
    localB->B_41_0_0.re = (rtu_In.re + rtu_In_f.re) + rtu_In_o.re;
    localB->B_41_0_0.im = (rtu_In.im + rtu_In_f.im) + rtu_In_o.im;
    localB->B_41_1_0.re = localP->P_0 * localB->B_41_0_0.re;
    localB->B_41_1_0.im = localP->P_0 * localB->B_41_0_0.im;
    if (ssIsMajorTimeStep(S) != 0) {
        srUpdateBC(localDW->ZeroSeqComputation_SubsysRanBC);
    }
}






      
  

        
    
      
          }
      
  



                





            /* 
 * Termination for enable system:
 *    '<S659>/Zero Seq. Computation'
 *    '<S701>/Zero Seq. Computation'
 *    '<S743>/Zero Seq. Computation'
 *    '<S785>/Zero Seq. Computation'
 */

          
  
  
      
  
      void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_ZeroSeqComputation_Term(SimStruct *const S)
  {
  
          
    
    
  

                
  

                  
  

          }
      
  



      
                                                                                                                                                                                                                                                                                
      
       
    
  

          /* Outputs for root system: '<Root>' */
      
            
  
  
      
  
     static  void mdlOutputs(SimStruct *S, int_T tid)
  {
  
                          /* local block i/o variables */
                    
                creal_T B_97_1340_0;


    
                creal_T B_97_1347_0;


    
                creal_T B_97_1354_0;


    
                creal_T B_97_1551_0;


    
                creal_T B_97_1558_0;


    
                creal_T B_97_1565_0;


    
                real_T B_97_6_0[5];


    
                real_T B_97_37_0[5];


    
                real_T B_97_80_0;


    
                real_T B_97_91_0;


    
                real_T B_97_92_0;


    
                real_T B_97_110_0;


    
                real_T B_97_114_0;


    
                real_T B_97_122_0;


    
                real_T B_97_134_0;


    
                real_T B_97_163_0;


    
                real_T B_97_374_0[6];


    
                real_T B_97_395_0;


    
                real_T B_97_397_0;


    
                real_T B_97_418_0;


    
                real_T B_97_427_0[2];


    
                real_T B_97_439_0;


    
                real_T B_97_478_0;


    
                real_T B_97_479_0;


    
                real_T B_97_490_0;


    
                real_T B_97_492_0;


    
                real_T B_97_507_0[3];


    
                real_T B_97_509_0[3];


    
                real_T B_97_772_0[6];


    
                real_T B_97_804_0;


    
                real_T B_97_806_0;


    
                real_T B_97_813_0[2];


    
                real_T B_97_825_0;


    
                real_T B_97_864_0;


    
                real_T B_97_865_0;


    
                real_T B_97_876_0;


    
                real_T B_97_878_0;


    
                real_T B_97_1082_0;


    
                real_T B_97_1083_0;


    
                real_T B_97_1094_0;


    
                real_T B_97_1096_0;


    
                real_T B_97_1098_0;


    
                real_T B_97_1100_0;


    
                real_T B_97_1109_0;


    
                real_T B_97_1118_0;


    
                real_T B_97_1197_0;


    
                real_T B_97_1200_0;


    
                real_T B_97_1202_0;


    
                real_T B_97_1228_0;


    
                real_T B_97_1230_0;


    
                real_T B_97_1249_0;


    
                real_T B_97_1335_0;


    
                real_T B_97_1338_0;


    
                real_T B_97_1342_0;


    
                real_T B_97_1345_0;


    
                real_T B_97_1349_0;


    
                real_T B_97_1352_0;


    
                real_T B_97_1376_0[2];


    
                real_T B_97_1388_0;


    
                real_T B_97_1427_0;


    
                real_T B_97_1428_0;


    
                real_T B_97_1439_0;


    
                real_T B_97_1441_0;


    
                real_T B_97_1546_0;


    
                real_T B_97_1549_0;


    
                real_T B_97_1553_0;


    
                real_T B_97_1556_0;


    
                real_T B_97_1560_0;


    
                real_T B_97_1563_0;


    
                real_T B_97_1587_0[2];


    
                real_T B_97_1599_0;


    
                real_T B_97_1637_0;


    
                real_T B_97_1638_0;


    
                real_T B_97_1649_0;


    
                real_T B_97_1651_0;


    
                real_T B_97_1714_0;


    
                real_T B_97_1764_0;


    
                real_T B_97_1767_0;


    
                real_T B_97_1769_0;


    
                real_T B_97_1801_0;


    
                real_T B_93_36_0;


    
                real_T B_93_38_0;


    
                real_T B_77_36_0;


    
                real_T B_77_38_0;


    
                real_T B_66_36_0;


    
                real_T B_66_38_0;


    
                real_T B_55_36_0;


    
                real_T B_55_38_0;


    
                real_T B_44_36_0;


    
                real_T B_44_38_0;


    
                real_T B_36_48_0;


    
                real_T B_36_49_0;


    
                real_T B_36_87_0;


    
                real_T B_36_88_0;


    
                real_T B_36_99_0;


    
                real_T B_36_101_0;


    
                real_T B_30_36_0;


    
                real_T B_30_38_0;


    
                real_T B_27_197_0;


    
                real_T B_27_198_0;


    
                real_T B_27_209_0;


    
                real_T B_27_211_0;


    
                real_T B_27_213_0;


    
                real_T B_27_215_0;


    
                real_T B_27_224_0;


    
                real_T B_27_233_0;


    
                real_T B_23_36_0;


    
                real_T B_23_38_0;


    
                boolean_T B_97_496_0[3];


    
                boolean_T B_97_882_0[3];



        
        
        
        

  
    
real_T *lastU;
real_T rtb_B_97_5_0;
real_T rtb_B_97_5_1;
real_T rtb_B_97_36_0;
real_T rtb_B_97_36_1;
real_T rtb_B_97_1_0;
real_T rtb_B_97_16_0;
real_T rtb_B_97_25_0[5];
real_T rtb_B_97_27_0;
real_T rtb_B_97_28_0;
real_T rtb_B_97_30_0[3];
real_T rtb_B_97_32_0;
real_T rtb_B_97_47_0;
real_T rtb_B_97_56_0[5];
real_T rtb_B_97_66_0;
real_T rtb_B_97_67_0;
real_T rtb_B_97_68_0;
boolean_T rtb_B_97_76_0;
real_T rtb_B_97_83_0;
real_T rtb_B_97_84_0;
real_T rtb_B_97_103_0;
real_T rtb_B_97_82_0;
boolean_T rtb_B_97_111_0;
real_T rtb_B_97_145_0;
uint8_T rtb_B_97_160_0;
real_T rtb_B_97_176_0;
real_T rtb_B_97_177_0;
real_T rtb_B_97_182_0;
real_T rtb_B_97_185_0;
real_T rtb_B_97_215_0;
real_T rtb_B_97_216_0;
real_T rtb_B_97_221_0;
real_T rtb_B_97_224_0;
real_T tmpForInput[9];
real_T rtb_B_97_396_0;
real_T rtb_B_97_409_0;
real_T rtb_B_97_410_0;
real_T rtb_B_97_411_0;
int8_T rtPrevAction;
int8_T rtAction;
real_T rtb_B_97_574_0;
real_T rtb_B_97_575_0;
real_T rtb_B_97_580_0;
real_T rtb_B_97_583_0;
real_T rtb_B_97_613_0;
real_T rtb_B_97_614_0;
real_T rtb_B_97_619_0;
real_T rtb_B_97_622_0;
real_T rtb_B_97_777_0;
real_T rtb_B_97_803_0;
real_T rtb_B_97_805_0;
real_T rtb_B_97_966_0;
real_T rtb_B_97_967_0;
real_T rtb_B_97_972_0;
real_T rtb_B_97_975_0;
real_T rtb_B_97_986_0;
real_T rtb_B_97_989_0;
real_T rtb_B_97_1173_0;
real_T rtb_B_97_1177_0[2];
real_T rtb_B_97_1194_0[18];
real_T rtb_B_97_17_0;
real_T rtb_B_97_1236_0;
real_T rtb_B_97_48_0;
real_T rtb_B_97_1718_0;
real_T rtb_B_97_1740_0[2];
real_T rtb_B_97_1804_0;
real_T rtb_B_27_80_0;
real_T rtb_B_27_81_0;
real_T rtb_B_27_86_0;
real_T rtb_B_27_89_0;
real_T rtb_B_27_100_0;
real_T rtb_B_27_103_0;
real_T rtb_B_77_37_0;
real_T rtb_B_97_1948_0;
real_T rtb_B_97_1997_0;
real_T rtb_B_97_2046_0;
real_T rtb_B_97_186_0;
real_T rtb_B_97_187_0;
real_T rtb_B_97_225_0;
real_T rtb_B_97_226_0;
real_T rtb_B_97_584_0;
real_T rtb_B_97_585_0;
real_T rtb_B_97_623_0;
real_T rtb_B_97_624_0;
real_T rtb_B_97_976_0;
real_T rtb_B_97_977_0;
real_T rtb_B_97_990_0;
real_T rtb_B_97_991_0;
real_T rtb_B_27_90_0;
real_T rtb_B_27_91_0;
real_T rtb_B_27_104_0;
real_T rtb_B_27_105_0;
real_T rtb_B_87_4_0[25];
int32_T i;
real_T tmp[25];
real_T tmp_0[5];
int32_T i_0;
boolean_T rtb_B_97_384_0;
real_T B_97_412_0_idx_0;
real_T B_97_412_0_idx_1;
real_T B_97_412_0_idx_2;
real_T B_5_0_0_idx_3;
real_T B_97_173_0_idx_0;
real_T B_97_173_0_idx_1;
real_T B_97_173_0_idx_2;
real_T B_97_212_0_idx_0;
real_T B_97_212_0_idx_1;
real_T B_97_212_0_idx_2;
real_T B_97_321_0_idx_0;
real_T B_97_321_0_idx_1;
real_T B_97_321_0_idx_2;
boolean_T B_97_384_0_idx_1;
boolean_T B_97_385_0_idx_1;
real_T B_97_328_1_idx_1;
real_T B_97_328_1_idx_2;
real_T B_97_719_0_idx_0;
real_T B_97_719_0_idx_1;
real_T B_97_719_0_idx_2;
real_T B_27_79_0_idx_0;
real_T B_27_79_0_idx_1;
real_T B_27_79_0_idx_2;
real_T B_27_95_0_idx_0;
real_T B_27_95_0_idx_1;
real_T B_27_95_0_idx_2;
real_T B_97_726_1_idx_0;
real_T B_97_726_1_idx_1;
real_T B_97_726_1_idx_2;
real_T B_97_965_0_idx_0;
real_T B_97_965_0_idx_1;
real_T B_97_965_0_idx_2;
real_T B_97_1542_0_idx_0;
real_T B_97_1542_0_idx_1;
real_T B_97_1542_0_idx_2;
real_T B_97_328_0_idx_0;
real_T B_97_328_0_idx_1;
real_T B_97_328_0_idx_2;
real_T B_97_291_1_idx_2;
real_T B_97_726_0_idx_0;
real_T B_97_726_0_idx_1;
real_T B_97_726_0_idx_2;
real_T B_97_781_0_idx_0;
real_T B_97_781_0_idx_1;
real_T B_97_781_0_idx_2;
real_T B_97_1331_0_idx_0;
real_T B_97_1331_0_idx_1;
real_T B_97_1331_0_idx_2;
real_T B_97_1542_1_idx_1;
real_T B_97_1542_1_idx_2;
B_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *_rtB;
P_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *_rtP;
X_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *_rtX;
DW_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *_rtDW;

    
  

              
  
    
  
    
  

    
                   
        
                  














































































































_rtDW = ((DW_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetRootDWork(S));
_rtX = ((X_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetContStates(S));
_rtP = ((P_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetModelRtp(S));
_rtB = ((B_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) _ssGetModelBlockIO(S));
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    rtb_B_97_1_0 = _rtDW->DiscreteTimeIntegrator1_DSTATE;
    muDoubleScalarSinCos(_rtDW->DiscreteTimeIntegrator1_DSTATE + _rtP->P_540 * ssGetTaskTime(S, 2), &rtb_B_97_5_0, &rtb_B_97_5_1);
    for (i = 0; i < 5; i++) {
        B_97_6_0[i] = _rtDW->fluxes_DSTATE[i];
    }
    _rtB->B_97_8_0 = _rtB->B_97_7_0;



      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_Saturation(S, _rtB->B_97_8_0, B_97_6_0, &_rtB->Saturation_i, &_rtDW->Saturation_i, &_rtP->Saturation_i);
    rtb_B_97_16_0 = _rtDW->DiscreteTimeIntegrator_DSTATE;
    rtb_B_97_17_0 = _rtB->B_97_15_0 + _rtDW->DiscreteTimeIntegrator_DSTATE;
    if (ssGetTaskTime(S, 2) >= _rtP->P_555) {
        memcpy(&rtb_B_87_4_0[0], &_rtB->B_97_14_0[0], 25U * sizeof(real_T));
        rtb_B_87_4_0[5] = rtb_B_97_17_0;
        rtb_B_87_4_0[1] = _rtP->P_498 * rtb_B_97_17_0;
        rtb_B_97_384_0 = (_rtB->B_97_18_0 >= _rtP->P_497);
        for (i = 0; i < 25; i++) {
            if (rtb_B_97_384_0) {
                B_5_0_0_idx_3 = _rtB->Saturation_i.B_3_28_0[i];
            } else {
                B_5_0_0_idx_3 = _rtB->B_97_19_0[i];
            }
            tmp[i] = ((0.0 - rtb_B_87_4_0[i]) - B_5_0_0_idx_3) * _rtP->P_499 + _rtB->B_97_13_0[i];
        }
        for (i = 0; i < 5; i++) {
            tmp_0[i] = 0.0;
            for (i_0 = 0; i_0 < 5; i_0++) {
                tmp_0[i] += tmp[5 * i_0 + i] * B_97_6_0[i_0];
            }
        }
        for (i = 0; i < 5; i++) {
            rtb_B_97_56_0[i] = 0.0;
            for (i_0 = 0; i_0 < 5; i_0++) {
                rtb_B_97_56_0[i] += _rtB->B_97_20_0[5 * i_0 + i] * _rtDW->voltages_DSTATE[i_0];
            }
            _rtB->B_97_24_0[i] = tmp_0[i] + rtb_B_97_56_0[i];
        }
    } else {
        for (i = 0; i < 5; i++) {
            _rtB->B_97_24_0[i] = B_97_6_0[i];
        }
    }
    rtb_B_97_384_0 = (_rtB->B_97_10_0 >= _rtP->P_545);
    for (i = 0; i < 25; i++) {
        if (rtb_B_97_384_0) {
            tmp[i] = _rtB->Saturation_i.B_3_27_0[i];
        } else {
            tmp[i] = _rtB->B_97_11_0[i];
        }
    }
    for (i = 0; i < 5; i++) {
        tmp_0[i] = 0.0;
        for (i_0 = 0; i_0 < 5; i_0++) {
            tmp_0[i] += tmp[5 * i_0 + i] * _rtB->B_97_24_0[i_0];
        }
        rtb_B_97_25_0[i] = _rtP->P_556[i] * tmp_0[i];
    }
    rtb_B_97_30_0[0] = _rtP->P_557 * (rtb_B_97_25_0[0] * rtb_B_97_5_1 + rtb_B_97_25_0[1] * rtb_B_97_5_0);
    rtb_B_97_30_0[1] = _rtP->P_557 * (((-rtb_B_97_25_0[0] - 1.7320508075688772 * rtb_B_97_25_0[1]) * rtb_B_97_5_1 + (1.7320508075688772 * rtb_B_97_25_0[0] - rtb_B_97_25_0[1]) * rtb_B_97_5_0) * 0.5);
    i = ssIsSampleHit(S, 2, 0);
    if (i != 0) {
        _rtB->B_97_31_0[0] = rtb_B_97_30_0[0];
        _rtB->B_97_31_0[1] = rtb_B_97_30_0[1];
    }
    rtb_B_97_32_0 = _rtDW->DiscreteTimeIntegrator1_DSTATE_h;
    muDoubleScalarSinCos(_rtDW->DiscreteTimeIntegrator1_DSTATE_h + _rtP->P_560 * ssGetTaskTime(S, 2), &rtb_B_97_36_0, &rtb_B_97_36_1);
    for (i = 0; i < 5; i++) {
        B_97_37_0[i] = _rtDW->fluxes_DSTATE_i[i];
    }
    _rtB->B_97_39_0 = _rtB->B_97_38_0;



      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_Saturation(S, _rtB->B_97_39_0, B_97_37_0, &_rtB->Saturation, &_rtDW->Saturation, &_rtP->Saturation);
    rtb_B_97_47_0 = _rtDW->DiscreteTimeIntegrator_DSTATE_f;
    rtb_B_97_48_0 = _rtB->B_97_46_0 + _rtDW->DiscreteTimeIntegrator_DSTATE_f;
    if (ssGetTaskTime(S, 2) >= _rtP->P_575) {
        memcpy(&rtb_B_87_4_0[0], &_rtB->B_97_45_0[0], 25U * sizeof(real_T));
        rtb_B_87_4_0[5] = rtb_B_97_48_0;
        rtb_B_87_4_0[1] = _rtP->P_1 * rtb_B_97_48_0;
        rtb_B_97_384_0 = (_rtB->B_97_49_0 >= _rtP->P_0);
        for (i = 0; i < 25; i++) {
            if (rtb_B_97_384_0) {
                B_5_0_0_idx_3 = _rtB->Saturation.B_3_28_0[i];
            } else {
                B_5_0_0_idx_3 = _rtB->B_97_50_0[i];
            }
            tmp[i] = ((0.0 - rtb_B_87_4_0[i]) - B_5_0_0_idx_3) * _rtP->P_2 + _rtB->B_97_44_0[i];
        }
        for (i = 0; i < 5; i++) {
            tmp_0[i] = 0.0;
            for (i_0 = 0; i_0 < 5; i_0++) {
                tmp_0[i] += tmp[5 * i_0 + i] * B_97_37_0[i_0];
            }
        }
        for (i = 0; i < 5; i++) {
            rtb_B_97_56_0[i] = 0.0;
            for (i_0 = 0; i_0 < 5; i_0++) {
                rtb_B_97_56_0[i] += _rtB->B_97_51_0[5 * i_0 + i] * _rtDW->voltages_DSTATE_c[i_0];
            }
            _rtB->B_97_55_0[i] = tmp_0[i] + rtb_B_97_56_0[i];
        }
    } else {
        for (i = 0; i < 5; i++) {
            _rtB->B_97_55_0[i] = B_97_37_0[i];
        }
    }
    rtb_B_97_384_0 = (_rtB->B_97_41_0 >= _rtP->P_565);
    for (i = 0; i < 25; i++) {
        if (rtb_B_97_384_0) {
            tmp[i] = _rtB->Saturation.B_3_27_0[i];
        } else {
            tmp[i] = _rtB->B_97_42_0[i];
        }
    }
    for (i = 0; i < 5; i++) {
        tmp_0[i] = 0.0;
        for (i_0 = 0; i_0 < 5; i_0++) {
            tmp_0[i] += tmp[5 * i_0 + i] * _rtB->B_97_55_0[i_0];
        }
        rtb_B_97_56_0[i] = _rtP->P_576[i] * tmp_0[i];
    }
    rtb_B_97_30_0[0] = _rtP->P_577 * (rtb_B_97_56_0[0] * rtb_B_97_36_1 + rtb_B_97_56_0[1] * rtb_B_97_36_0);
    rtb_B_97_30_0[1] = _rtP->P_577 * (((-rtb_B_97_56_0[0] - 1.7320508075688772 * rtb_B_97_56_0[1]) * rtb_B_97_36_1 + (1.7320508075688772 * rtb_B_97_56_0[0] - rtb_B_97_56_0[1]) * rtb_B_97_36_0) * 0.5);
    i = ssIsSampleHit(S, 2, 0);
    if (i != 0) {
        _rtB->B_97_62_0[0] = rtb_B_97_30_0[0];
        _rtB->B_97_62_0[1] = rtb_B_97_30_0[1];
    }
}
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    rtb_B_97_66_0 = _rtDW->itinit1_PreviousInput;
    rtb_B_97_67_0 = _rtP->P_582 * _rtDW->itinit1_PreviousInput;
    rtb_B_97_68_0 = 1.000001 * rtb_B_97_67_0 * 0.96711798839458663 / 0.9999;
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_69_0 = _rtP->P_583 * _rtDW->Currentfilter_states;
    _rtB->B_97_72_0 = (_rtB->B_97_69_0 > Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_rtC(S)->B_97_70_0);
}
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    _rtB->B_97_73_0 = _rtDW->itinit_PreviousInput;
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    if (_rtDW->inti_IC_LOADING != 0) {
        _rtDW->inti_DSTATE = _rtB->B_97_73_0;
        if (_rtDW->inti_DSTATE >= _rtP->P_589) {
            _rtDW->inti_DSTATE = _rtP->P_589;
        } else {
            if (_rtDW->inti_DSTATE <= _rtP->P_590) {
                _rtDW->inti_DSTATE = _rtP->P_590;
            }
        }
    }
    if ((_rtB->B_97_72_0 > 0.0) && (_rtDW->inti_PrevResetState <= 0)) {
        _rtDW->inti_DSTATE = _rtB->B_97_73_0;
        if (_rtDW->inti_DSTATE >= _rtP->P_589) {
            _rtDW->inti_DSTATE = _rtP->P_589;
        } else {
            if (_rtDW->inti_DSTATE <= _rtP->P_590) {
                _rtDW->inti_DSTATE = _rtP->P_590;
            }
        }
    }
    if (_rtDW->inti_DSTATE >= _rtP->P_589) {
        _rtDW->inti_DSTATE = _rtP->P_589;
    } else {
        if (_rtDW->inti_DSTATE <= _rtP->P_590) {
            _rtDW->inti_DSTATE = _rtP->P_590;
        }
    }
    _rtB->B_97_75_0 = _rtP->P_591 * _rtDW->inti_DSTATE;
}
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    rtb_B_97_76_0 = (_rtB->B_97_75_0 > rtb_B_97_67_0);
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    if (_rtB->B_97_75_0 < _rtB->B_97_77_0) {
        _rtB->B_97_79_0 = _rtB->B_97_77_0;
    } else {
        _rtB->B_97_79_0 = _rtB->B_97_75_0;
    }
}
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    if (rtb_B_97_76_0) {
        B_97_80_0 = rtb_B_97_67_0;
    } else {
        B_97_80_0 = _rtB->B_97_79_0;
    }
    if ((rtb_B_97_68_0 <= B_97_80_0) > _rtP->P_593) {
        rtb_B_97_82_0 = rtb_B_97_67_0;
    } else {
        rtb_B_97_82_0 = B_97_80_0;
    }
    rtb_B_97_83_0 = -0.027082910519425695 * rtb_B_97_66_0 / (rtb_B_97_66_0 - rtb_B_97_82_0) * rtb_B_97_82_0;
    rtb_B_97_84_0 = -_rtB->B_97_72_0 * 0.027082910519425695 * _rtB->B_97_69_0 * rtb_B_97_66_0 / (rtb_B_97_66_0 - rtb_B_97_82_0);
    rtb_B_97_67_0 = _rtP->P_596 * rtb_B_97_66_0;
    rtb_B_97_68_0 = -rtb_B_97_67_0 * 0.999 * 0.1 * 0.9999;
    if (_rtB->B_97_75_0 < rtb_B_97_68_0) {
        B_97_91_0 = rtb_B_97_68_0;
    } else {
        B_97_91_0 = _rtB->B_97_75_0;
    }
    if (_rtB->B_97_75_0 > rtb_B_97_67_0) {
        B_97_92_0 = rtb_B_97_67_0;
    } else {
        B_97_92_0 = B_97_91_0;
    }
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_95_0 = (_rtB->B_97_69_0 < Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_rtC(S)->B_97_93_0);
}
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    switch ((int32_T)_rtB->B_97_85_0) {
      case 1:
        rtb_B_97_103_0 = -(_rtB->B_97_86_0 * _rtB->B_97_95_0) * 0.027082910519425695 * (_rtB->B_97_86_0 * _rtB->B_97_69_0) * (206.79999999999947 / (_rtB->B_97_86_0 * B_97_92_0 + 20.67999999999995));
        break;
      case 2:
        B_5_0_0_idx_3 = _rtB->B_97_96_0 * rtb_B_97_66_0;
        rtb_B_97_103_0 = -(_rtB->B_97_96_0 * _rtB->B_97_95_0) * 0.027082910519425695 * (_rtB->B_97_96_0 * _rtB->B_97_69_0) * B_5_0_0_idx_3 / (_rtB->B_97_96_0 * B_97_92_0 + B_5_0_0_idx_3 * 0.1);
        break;
      case 3:
        rtb_B_97_103_0 = -(_rtB->B_97_97_0 * _rtB->B_97_95_0) * 0.027082910519425695 * (_rtB->B_97_97_0 * _rtB->B_97_69_0) * (206.79999999999947 / (muDoubleScalarAbs(_rtB->B_97_97_0 * B_97_92_0) + 20.67999999999995));
        break;
      default:
        rtb_B_97_103_0 = -(_rtB->B_97_98_0 * _rtB->B_97_95_0) * 0.027082910519425695 * (_rtB->B_97_98_0 * _rtB->B_97_69_0) * (206.79999999999947 / (muDoubleScalarAbs(_rtB->B_97_98_0 * B_97_92_0) + 20.67999999999995));
        break;
    }
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_105_0 = _rtDW->DiscreteTimeIntegrator_DSTATE_l;
}
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    switch ((int32_T)_rtB->B_97_104_0) {
      case 1:
        rtb_B_97_67_0 = _rtB->B_97_105_0;
        break;
      case 2:
        if (rtb_B_97_82_0 > _rtP->P_3) {
            B_5_0_0_idx_3 = _rtP->P_3;
        } else if (rtb_B_97_82_0 < _rtP->P_4) {
            B_5_0_0_idx_3 = _rtP->P_4;
        } else {
            B_5_0_0_idx_3 = rtb_B_97_82_0;
        }
        rtb_B_97_67_0 = muDoubleScalarExp(-0.30530973451327431 * B_5_0_0_idx_3) * 60.7146212936532;
        break;
      case 3:
        rtb_B_97_67_0 = _rtB->B_97_105_0;
        break;
      default:
        rtb_B_97_67_0 = _rtB->B_97_105_0;
        break;
    }
    B_97_110_0 = ((((rtb_B_97_83_0 + rtb_B_97_84_0) + rtb_B_97_103_0) + rtb_B_97_67_0) + -0.0 * rtb_B_97_82_0) + _rtB->B_97_65_0;
    rtb_B_97_111_0 = (B_97_110_0 > _rtB->B_97_64_0);
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_112_0 = _rtDW->Memory2_PreviousInput;
}
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    if (B_97_110_0 < _rtB->B_97_112_0) {
        B_97_114_0 = _rtB->B_97_112_0;
    } else {
        B_97_114_0 = B_97_110_0;
    }
    if (rtb_B_97_111_0) {
        _rtB->B_97_115_0 = _rtB->B_97_64_0;
    } else {
        _rtB->B_97_115_0 = B_97_114_0;
    }
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
                
          /* Level2 S-Function Block: '<S1410>/B_82_9' (sfun_spssw_discc) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 116, SS_CALL_MDL_OUTPUTS);

   
  
    B_97_122_0 = _rtP->P_615 * _rtB->B_97_116_0[57] * (_rtP->P_616 * _rtB->B_97_116_0[38]);
    _rtB->B_97_129_0 = _rtP->P_617 * _rtB->B_97_116_0[58];
    _rtB->B_97_133_0 = _rtP->P_618 * _rtB->B_97_129_0;
}
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    B_97_134_0 = _rtB->B_97_115_0 - _rtB->B_97_133_0;
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_143_0 = ((real_T)(_rtB->B_97_129_0 < Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_rtC(S)->B_97_137_0) * _rtP->P_621 - _rtB->B_97_105_0) * muDoubleScalarAbs(_rtB->B_97_129_0) * _rtP->P_622;
}
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    rtb_B_97_66_0 *= _rtP->P_623;
    rtb_B_97_145_0 = (1.0 - rtb_B_97_82_0 / rtb_B_97_66_0) * 100.0;
    _rtB->B_97_146_0 = _rtP->P_624 * rtb_B_97_82_0;
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_147_0 = _rtP->P_625 * _rtB->B_97_129_0;
}
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    if (rtb_B_97_145_0 > _rtP->P_626) {
        rtb_B_97_145_0 = _rtP->P_626;
    } else {
        if (rtb_B_97_145_0 < _rtP->P_627) {
            rtb_B_97_145_0 = _rtP->P_627;
        }
    }
    rtb_B_97_160_0 = (uint8_T)(rtb_B_97_145_0 >= _rtB->B_97_159_0);
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    if (_rtB->B_97_116_1[0] >= _rtP->P_631) {
        B_97_163_0 = _rtP->P_5 * _rtB->B_97_116_0[0];
    } else {
        B_97_163_0 = _rtB->B_97_161_0;
    }
}
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    _rtB->B_97_164_0 = rtb_B_97_160_0;
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    B_97_173_0_idx_0 = _rtP->P_632 * _rtB->B_97_116_0[32] * _rtP->P_635;
    B_97_173_0_idx_1 = _rtP->P_633 * _rtB->B_97_116_0[33] * _rtP->P_635;
    B_97_173_0_idx_2 = _rtP->P_634 * _rtB->B_97_116_0[34] * _rtP->P_635;
    rtb_B_97_82_0 = _rtP->P_636 * ssGetTaskTime(S, 2);
    rtb_B_97_176_0 = muDoubleScalarSin(rtb_B_97_82_0);
    rtb_B_97_177_0 = muDoubleScalarCos(rtb_B_97_82_0);
    rtb_B_97_182_0 = (0.0 - rtb_B_97_176_0 * _rtB->B_97_178_0) - rtb_B_97_177_0 * _rtB->B_97_180_0;
    rtb_B_97_185_0 = rtb_B_97_176_0 * _rtB->B_97_180_0 - rtb_B_97_177_0 * _rtB->B_97_178_0;
    rtb_B_97_186_0 = (0.0 - rtb_B_97_182_0) - rtb_B_97_176_0;
    rtb_B_97_187_0 = (0.0 - rtb_B_97_185_0) - rtb_B_97_177_0;
    _rtB->B_97_188_0 = ((B_97_173_0_idx_0 * rtb_B_97_176_0 + B_97_173_0_idx_1 * rtb_B_97_182_0) + B_97_173_0_idx_2 * rtb_B_97_186_0) * 0.66666666666666663;
    if (_rtDW->Integ4_SYSTEM_ENABLE != 0) {
        _rtB->B_97_189_0 = _rtDW->Integ4_DSTATE;
    } else {
        _rtB->B_97_189_0 = _rtP->P_639 * _rtB->B_97_188_0 + _rtDW->Integ4_DSTATE;
    }
    _rtB->B_97_190_0 = _rtP->P_641;
                
          /* Level2 S-Function Block: '<S168>/B_82_31' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 191, SS_CALL_MDL_OUTPUTS);

   
  
    rtb_B_97_82_0 = ssGetTaskTime(S, 2);
    if (rtb_B_97_82_0 < _rtP->P_644) {
        rtb_B_97_82_0 = _rtP->P_645;
    } else {
        rtb_B_97_82_0 = _rtP->P_646;
    }
    if (rtb_B_97_82_0 >= _rtP->P_648) {
        rtb_B_97_66_0 = (_rtB->B_97_189_0 - _rtB->B_97_191_0) * _rtP->P_642 + (_rtP->P_11 * _rtB->B_97_188_0 - _rtP->P_10 * _rtDW->UnitDelay_DSTATE);
    } else {
        rtb_B_97_66_0 = _rtP->P_647;
    }
    _rtB->B_97_198_0 = ((B_97_173_0_idx_0 * rtb_B_97_177_0 + B_97_173_0_idx_1 * rtb_B_97_185_0) + B_97_173_0_idx_2 * rtb_B_97_187_0) * 0.66666666666666663;
    if (_rtDW->Integ4_SYSTEM_ENABLE_a != 0) {
        _rtB->B_97_199_0 = _rtDW->Integ4_DSTATE_f;
    } else {
        _rtB->B_97_199_0 = _rtP->P_649 * _rtB->B_97_198_0 + _rtDW->Integ4_DSTATE_f;
    }
    _rtB->B_97_200_0 = _rtP->P_651;
                
          /* Level2 S-Function Block: '<S169>/B_82_33' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 201, SS_CALL_MDL_OUTPUTS);

   
  
    rtb_B_97_82_0 = ssGetTaskTime(S, 2);
    if (rtb_B_97_82_0 < _rtP->P_654) {
        rtb_B_97_82_0 = _rtP->P_655;
    } else {
        rtb_B_97_82_0 = _rtP->P_656;
    }
    if (rtb_B_97_82_0 >= _rtP->P_658) {
        rtb_B_97_83_0 = (_rtB->B_97_199_0 - _rtB->B_97_201_0) * _rtP->P_652 + (_rtP->P_13 * _rtB->B_97_198_0 - _rtP->P_12 * _rtDW->UnitDelay_DSTATE_d);
    } else {
        rtb_B_97_83_0 = _rtP->P_657;
    }
    B_97_212_0_idx_0 = _rtP->P_659 * _rtB->B_97_116_0[54] * _rtP->P_662;
    B_97_212_0_idx_1 = _rtP->P_660 * _rtB->B_97_116_0[55] * _rtP->P_662;
    B_97_212_0_idx_2 = _rtP->P_661 * _rtB->B_97_116_0[56] * _rtP->P_662;
    rtb_B_97_82_0 = _rtP->P_663 * ssGetTaskTime(S, 2);
    rtb_B_97_215_0 = muDoubleScalarSin(rtb_B_97_82_0);
    rtb_B_97_216_0 = muDoubleScalarCos(rtb_B_97_82_0);
    rtb_B_97_221_0 = (0.0 - rtb_B_97_215_0 * _rtB->B_97_217_0) - rtb_B_97_216_0 * _rtB->B_97_219_0;
    rtb_B_97_224_0 = rtb_B_97_215_0 * _rtB->B_97_219_0 - rtb_B_97_216_0 * _rtB->B_97_217_0;
    rtb_B_97_225_0 = (0.0 - rtb_B_97_221_0) - rtb_B_97_215_0;
    rtb_B_97_226_0 = (0.0 - rtb_B_97_224_0) - rtb_B_97_216_0;
    _rtB->B_97_227_0 = ((B_97_212_0_idx_0 * rtb_B_97_215_0 + B_97_212_0_idx_1 * rtb_B_97_221_0) + B_97_212_0_idx_2 * rtb_B_97_225_0) * 0.66666666666666663;
    if (_rtDW->Integ4_SYSTEM_ENABLE_c != 0) {
        _rtB->B_97_228_0 = _rtDW->Integ4_DSTATE_m;
    } else {
        _rtB->B_97_228_0 = _rtP->P_666 * _rtB->B_97_227_0 + _rtDW->Integ4_DSTATE_m;
    }
    _rtB->B_97_229_0 = _rtP->P_668;
                
          /* Level2 S-Function Block: '<S162>/B_82_35' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 230, SS_CALL_MDL_OUTPUTS);

   
  
    rtb_B_97_82_0 = ssGetTaskTime(S, 2);
    if (rtb_B_97_82_0 < _rtP->P_671) {
        rtb_B_97_82_0 = _rtP->P_672;
    } else {
        rtb_B_97_82_0 = _rtP->P_673;
    }
    if (rtb_B_97_82_0 >= _rtP->P_675) {
        rtb_B_97_84_0 = (_rtB->B_97_228_0 - _rtB->B_97_230_0) * _rtP->P_669 + (_rtP->P_7 * _rtB->B_97_227_0 - _rtP->P_6 * _rtDW->UnitDelay_DSTATE_e);
    } else {
        rtb_B_97_84_0 = _rtP->P_674;
    }
    _rtB->B_97_237_0 = ((B_97_212_0_idx_0 * rtb_B_97_216_0 + B_97_212_0_idx_1 * rtb_B_97_224_0) + B_97_212_0_idx_2 * rtb_B_97_226_0) * 0.66666666666666663;
    if (_rtDW->Integ4_SYSTEM_ENABLE_l != 0) {
        _rtB->B_97_238_0 = _rtDW->Integ4_DSTATE_d;
    } else {
        _rtB->B_97_238_0 = _rtP->P_676 * _rtB->B_97_237_0 + _rtDW->Integ4_DSTATE_d;
    }
    _rtB->B_97_239_0 = _rtP->P_678;
                
          /* Level2 S-Function Block: '<S163>/B_82_37' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 240, SS_CALL_MDL_OUTPUTS);

   
  
    rtb_B_97_82_0 = ssGetTaskTime(S, 2);
    if (rtb_B_97_82_0 < _rtP->P_681) {
        rtb_B_97_82_0 = _rtP->P_682;
    } else {
        rtb_B_97_82_0 = _rtP->P_683;
    }
    if (rtb_B_97_82_0 >= _rtP->P_685) {
        rtb_B_97_103_0 = (_rtB->B_97_238_0 - _rtB->B_97_240_0) * _rtP->P_679 + (_rtP->P_9 * _rtB->B_97_237_0 - _rtP->P_8 * _rtDW->UnitDelay_DSTATE_f);
    } else {
        rtb_B_97_103_0 = _rtP->P_684;
    }
    rtb_B_97_82_0 = muDoubleScalarHypot(rtb_B_97_66_0, rtb_B_97_83_0) * muDoubleScalarHypot(rtb_B_97_84_0, rtb_B_97_103_0);
    rtb_B_97_66_0 = (_rtP->P_686 * muDoubleScalarAtan2(rtb_B_97_83_0, rtb_B_97_66_0) - _rtP->P_687 * muDoubleScalarAtan2(rtb_B_97_103_0, rtb_B_97_84_0)) * _rtP->P_688;
    _rtB->B_97_259_0[0] = rtb_B_97_82_0 * muDoubleScalarCos(rtb_B_97_66_0) * _rtP->P_689;
    _rtB->B_97_259_0[1] = muDoubleScalarSin(rtb_B_97_66_0) * rtb_B_97_82_0 * _rtP->P_689;
    _rtB->B_97_260_0[0] = _rtP->P_690 * _rtB->B_97_259_0[0];
    _rtB->B_97_260_0[1] = _rtP->P_690 * _rtB->B_97_259_0[1];
    ssCallAccelRunBlock(S, 97, 261, SS_CALL_MDL_OUTPUTS);
    _rtB->B_97_265_0[0] = _rtP->P_691 * _rtB->B_97_116_0[48] * _rtP->P_694;
    _rtB->B_97_265_0[1] = _rtP->P_692 * _rtB->B_97_116_0[49] * _rtP->P_694;
    _rtB->B_97_265_0[2] = _rtP->P_693 * _rtB->B_97_116_0[50] * _rtP->P_694;
    if (_rtDW->systemEnable != 0) {
        _rtDW->lastSin = muDoubleScalarSin(_rtP->P_697 * ssGetTaskTime(S, 2));
        _rtDW->lastCos = muDoubleScalarCos(_rtP->P_697 * ssGetTaskTime(S, 2));
        _rtDW->systemEnable = 0;
    }
    rtb_B_97_82_0 = ((_rtDW->lastSin * _rtP->P_701 + _rtDW->lastCos * _rtP->P_700) * _rtP->P_699 + (_rtDW->lastCos * _rtP->P_701 - _rtDW->lastSin * _rtP->P_700) * _rtP->P_698) * _rtP->P_695 + _rtP->P_696;
    _rtB->B_97_267_0[0] = _rtB->B_97_265_0[0] * rtb_B_97_82_0;
    _rtB->B_97_267_0[1] = _rtB->B_97_265_0[1] * rtb_B_97_82_0;
    _rtB->B_97_267_0[2] = _rtB->B_97_265_0[2] * rtb_B_97_82_0;
    if (_rtDW->Integ4_SYSTEM_ENABLE_p != 0) {
        _rtB->B_97_268_0[0] = _rtDW->Integ4_DSTATE_k[0];
        _rtB->B_97_268_0[1] = _rtDW->Integ4_DSTATE_k[1];
        _rtB->B_97_268_0[2] = _rtDW->Integ4_DSTATE_k[2];
    } else {
        _rtB->B_97_268_0[0] = _rtP->P_702 * _rtB->B_97_267_0[0] + _rtDW->Integ4_DSTATE_k[0];
        _rtB->B_97_268_0[1] = _rtP->P_702 * _rtB->B_97_267_0[1] + _rtDW->Integ4_DSTATE_k[1];
        _rtB->B_97_268_0[2] = _rtP->P_702 * _rtB->B_97_267_0[2] + _rtDW->Integ4_DSTATE_k[2];
    }
    _rtB->B_97_269_0 = _rtP->P_704;
                
          /* Level2 S-Function Block: '<S694>/B_82_39' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 270, SS_CALL_MDL_OUTPUTS);

   
  
    if (ssGetTaskTime(S, 2) >= _rtB->B_97_269_0) {
        _rtB->B_97_277_0[0] = (_rtB->B_97_268_0[0] - _rtB->B_97_270_0[0]) * _rtP->P_713 + (_rtP->P_372 * _rtB->B_97_267_0[0] - _rtP->P_371 * _rtDW->UnitDelay_DSTATE_dz[0]);
        _rtB->B_97_277_0[1] = (_rtB->B_97_268_0[1] - _rtB->B_97_270_0[1]) * _rtP->P_713 + (_rtP->P_372 * _rtB->B_97_267_0[1] - _rtP->P_371 * _rtDW->UnitDelay_DSTATE_dz[1]);
        _rtB->B_97_277_0[2] = (_rtB->B_97_268_0[2] - _rtB->B_97_270_0[2]) * _rtP->P_713 + (_rtP->P_372 * _rtB->B_97_267_0[2] - _rtP->P_371 * _rtDW->UnitDelay_DSTATE_dz[2]);
    } else {
        _rtB->B_97_277_0[0] = _rtDW->UnitDelay1_DSTATE[0];
        _rtB->B_97_277_0[1] = _rtDW->UnitDelay1_DSTATE[1];
        _rtB->B_97_277_0[2] = _rtDW->UnitDelay1_DSTATE[2];
    }
    if (_rtDW->systemEnable_a != 0) {
        _rtDW->lastSin_m = muDoubleScalarSin(_rtP->P_718 * ssGetTaskTime(S, 2));
        _rtDW->lastCos_i = muDoubleScalarCos(_rtP->P_718 * ssGetTaskTime(S, 2));
        _rtDW->systemEnable_a = 0;
    }
    rtb_B_97_82_0 = ((_rtDW->lastSin_m * _rtP->P_722 + _rtDW->lastCos_i * _rtP->P_721) * _rtP->P_720 + (_rtDW->lastCos_i * _rtP->P_722 - _rtDW->lastSin_m * _rtP->P_721) * _rtP->P_719) * _rtP->P_716 + _rtP->P_717;
    _rtB->B_97_279_0[0] = _rtB->B_97_265_0[0] * rtb_B_97_82_0;
    _rtB->B_97_279_0[1] = _rtB->B_97_265_0[1] * rtb_B_97_82_0;
    _rtB->B_97_279_0[2] = _rtB->B_97_265_0[2] * rtb_B_97_82_0;
    if (_rtDW->Integ4_SYSTEM_ENABLE_f != 0) {
        _rtB->B_97_280_0[0] = _rtDW->Integ4_DSTATE_h[0];
        _rtB->B_97_280_0[1] = _rtDW->Integ4_DSTATE_h[1];
        _rtB->B_97_280_0[2] = _rtDW->Integ4_DSTATE_h[2];
    } else {
        _rtB->B_97_280_0[0] = _rtP->P_723 * _rtB->B_97_279_0[0] + _rtDW->Integ4_DSTATE_h[0];
        _rtB->B_97_280_0[1] = _rtP->P_723 * _rtB->B_97_279_0[1] + _rtDW->Integ4_DSTATE_h[1];
        _rtB->B_97_280_0[2] = _rtP->P_723 * _rtB->B_97_279_0[2] + _rtDW->Integ4_DSTATE_h[2];
    }
    _rtB->B_97_281_0 = _rtP->P_725;
                
          /* Level2 S-Function Block: '<S692>/B_82_41' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 282, SS_CALL_MDL_OUTPUTS);

   
  
    if (ssGetTaskTime(S, 2) >= _rtB->B_97_281_0) {
        _rtB->B_97_289_0[0] = (_rtB->B_97_280_0[0] - _rtB->B_97_282_0[0]) * _rtP->P_734 + (_rtP->P_370 * _rtB->B_97_279_0[0] - _rtP->P_369 * _rtDW->UnitDelay_DSTATE_o[0]);
        _rtB->B_97_289_0[1] = (_rtB->B_97_280_0[1] - _rtB->B_97_282_0[1]) * _rtP->P_734 + (_rtP->P_370 * _rtB->B_97_279_0[1] - _rtP->P_369 * _rtDW->UnitDelay_DSTATE_o[1]);
        _rtB->B_97_289_0[2] = (_rtB->B_97_280_0[2] - _rtB->B_97_282_0[2]) * _rtP->P_734 + (_rtP->P_370 * _rtB->B_97_279_0[2] - _rtP->P_369 * _rtDW->UnitDelay_DSTATE_o[2]);
    } else {
        _rtB->B_97_289_0[0] = _rtDW->UnitDelay1_DSTATE_c[0];
        _rtB->B_97_289_0[1] = _rtDW->UnitDelay1_DSTATE_c[1];
        _rtB->B_97_289_0[2] = _rtDW->UnitDelay1_DSTATE_c[2];
    }
    rtb_B_97_30_0[0] = muDoubleScalarHypot(_rtB->B_97_277_0[0], _rtB->B_97_289_0[0]);
    rtb_B_97_66_0 = muDoubleScalarAtan2(_rtB->B_97_289_0[0], _rtB->B_97_277_0[0]);
    rtb_B_97_30_0[1] = muDoubleScalarHypot(_rtB->B_97_277_0[1], _rtB->B_97_289_0[1]);
    B_5_0_0_idx_3 = muDoubleScalarAtan2(_rtB->B_97_289_0[1], _rtB->B_97_277_0[1]);
    rtb_B_97_30_0[2] = muDoubleScalarHypot(_rtB->B_97_277_0[2], _rtB->B_97_289_0[2]);
    B_97_291_1_idx_2 = muDoubleScalarAtan2(_rtB->B_97_289_0[2], _rtB->B_97_277_0[2]);
    _rtB->B_97_295_0[0] = _rtP->P_737 * _rtB->B_97_116_0[68] * _rtP->P_740;
    _rtB->B_97_295_0[1] = _rtP->P_738 * _rtB->B_97_116_0[69] * _rtP->P_740;
    _rtB->B_97_295_0[2] = _rtP->P_739 * _rtB->B_97_116_0[70] * _rtP->P_740;
    if (_rtDW->systemEnable_m != 0) {
        _rtDW->lastSin_e = muDoubleScalarSin(_rtP->P_743 * ssGetTaskTime(S, 2));
        _rtDW->lastCos_j = muDoubleScalarCos(_rtP->P_743 * ssGetTaskTime(S, 2));
        _rtDW->systemEnable_m = 0;
    }
    rtb_B_97_82_0 = ((_rtDW->lastSin_e * _rtP->P_747 + _rtDW->lastCos_j * _rtP->P_746) * _rtP->P_745 + (_rtDW->lastCos_j * _rtP->P_747 - _rtDW->lastSin_e * _rtP->P_746) * _rtP->P_744) * _rtP->P_741 + _rtP->P_742;
    _rtB->B_97_297_0[0] = _rtB->B_97_295_0[0] * rtb_B_97_82_0;
    _rtB->B_97_297_0[1] = _rtB->B_97_295_0[1] * rtb_B_97_82_0;
    _rtB->B_97_297_0[2] = _rtB->B_97_295_0[2] * rtb_B_97_82_0;
    if (_rtDW->Integ4_SYSTEM_ENABLE_h != 0) {
        _rtB->B_97_298_0[0] = _rtDW->Integ4_DSTATE_n[0];
        _rtB->B_97_298_0[1] = _rtDW->Integ4_DSTATE_n[1];
        _rtB->B_97_298_0[2] = _rtDW->Integ4_DSTATE_n[2];
    } else {
        _rtB->B_97_298_0[0] = _rtP->P_748 * _rtB->B_97_297_0[0] + _rtDW->Integ4_DSTATE_n[0];
        _rtB->B_97_298_0[1] = _rtP->P_748 * _rtB->B_97_297_0[1] + _rtDW->Integ4_DSTATE_n[1];
        _rtB->B_97_298_0[2] = _rtP->P_748 * _rtB->B_97_297_0[2] + _rtDW->Integ4_DSTATE_n[2];
    }
    _rtB->B_97_299_0 = _rtP->P_750;
                
          /* Level2 S-Function Block: '<S700>/B_82_43' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 300, SS_CALL_MDL_OUTPUTS);

   
  
    if (ssGetTaskTime(S, 2) >= _rtB->B_97_299_0) {
        _rtB->B_97_307_0[0] = (_rtB->B_97_298_0[0] - _rtB->B_97_300_0[0]) * _rtP->P_759 + (_rtP->P_376 * _rtB->B_97_297_0[0] - _rtP->P_375 * _rtDW->UnitDelay_DSTATE_ee[0]);
        _rtB->B_97_307_0[1] = (_rtB->B_97_298_0[1] - _rtB->B_97_300_0[1]) * _rtP->P_759 + (_rtP->P_376 * _rtB->B_97_297_0[1] - _rtP->P_375 * _rtDW->UnitDelay_DSTATE_ee[1]);
        _rtB->B_97_307_0[2] = (_rtB->B_97_298_0[2] - _rtB->B_97_300_0[2]) * _rtP->P_759 + (_rtP->P_376 * _rtB->B_97_297_0[2] - _rtP->P_375 * _rtDW->UnitDelay_DSTATE_ee[2]);
    } else {
        _rtB->B_97_307_0[0] = _rtDW->UnitDelay1_DSTATE_h[0];
        _rtB->B_97_307_0[1] = _rtDW->UnitDelay1_DSTATE_h[1];
        _rtB->B_97_307_0[2] = _rtDW->UnitDelay1_DSTATE_h[2];
    }
    if (_rtDW->systemEnable_l != 0) {
        _rtDW->lastSin_mo = muDoubleScalarSin(_rtP->P_764 * ssGetTaskTime(S, 2));
        _rtDW->lastCos_f = muDoubleScalarCos(_rtP->P_764 * ssGetTaskTime(S, 2));
        _rtDW->systemEnable_l = 0;
    }
    rtb_B_97_82_0 = ((_rtDW->lastSin_mo * _rtP->P_768 + _rtDW->lastCos_f * _rtP->P_767) * _rtP->P_766 + (_rtDW->lastCos_f * _rtP->P_768 - _rtDW->lastSin_mo * _rtP->P_767) * _rtP->P_765) * _rtP->P_762 + _rtP->P_763;
    _rtB->B_97_309_0[0] = _rtB->B_97_295_0[0] * rtb_B_97_82_0;
    _rtB->B_97_309_0[1] = _rtB->B_97_295_0[1] * rtb_B_97_82_0;
    _rtB->B_97_309_0[2] = _rtB->B_97_295_0[2] * rtb_B_97_82_0;
    if (_rtDW->Integ4_SYSTEM_ENABLE_k != 0) {
        _rtB->B_97_310_0[0] = _rtDW->Integ4_DSTATE_l[0];
        _rtB->B_97_310_0[1] = _rtDW->Integ4_DSTATE_l[1];
        _rtB->B_97_310_0[2] = _rtDW->Integ4_DSTATE_l[2];
    } else {
        _rtB->B_97_310_0[0] = _rtP->P_769 * _rtB->B_97_309_0[0] + _rtDW->Integ4_DSTATE_l[0];
        _rtB->B_97_310_0[1] = _rtP->P_769 * _rtB->B_97_309_0[1] + _rtDW->Integ4_DSTATE_l[1];
        _rtB->B_97_310_0[2] = _rtP->P_769 * _rtB->B_97_309_0[2] + _rtDW->Integ4_DSTATE_l[2];
    }
    _rtB->B_97_311_0 = _rtP->P_771;
                
          /* Level2 S-Function Block: '<S698>/B_82_45' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 312, SS_CALL_MDL_OUTPUTS);

   
  
    if (ssGetTaskTime(S, 2) >= _rtB->B_97_311_0) {
        _rtB->B_97_319_0[0] = (_rtB->B_97_310_0[0] - _rtB->B_97_312_0[0]) * _rtP->P_780 + (_rtP->P_374 * _rtB->B_97_309_0[0] - _rtP->P_373 * _rtDW->UnitDelay_DSTATE_p[0]);
        _rtB->B_97_319_0[1] = (_rtB->B_97_310_0[1] - _rtB->B_97_312_0[1]) * _rtP->P_780 + (_rtP->P_374 * _rtB->B_97_309_0[1] - _rtP->P_373 * _rtDW->UnitDelay_DSTATE_p[1]);
        _rtB->B_97_319_0[2] = (_rtB->B_97_310_0[2] - _rtB->B_97_312_0[2]) * _rtP->P_780 + (_rtP->P_374 * _rtB->B_97_309_0[2] - _rtP->P_373 * _rtDW->UnitDelay_DSTATE_p[2]);
    } else {
        _rtB->B_97_319_0[0] = _rtDW->UnitDelay1_DSTATE_j[0];
        _rtB->B_97_319_0[1] = _rtDW->UnitDelay1_DSTATE_j[1];
        _rtB->B_97_319_0[2] = _rtDW->UnitDelay1_DSTATE_j[2];
    }
    B_97_321_0_idx_0 = rtb_B_97_30_0[0] * muDoubleScalarHypot(_rtB->B_97_307_0[0], _rtB->B_97_319_0[0]) * _rtP->P_783;
    rtb_B_97_66_0 = (_rtP->P_784 * rtb_B_97_66_0 - _rtP->P_785 * muDoubleScalarAtan2(_rtB->B_97_319_0[0], _rtB->B_97_307_0[0])) * _rtP->P_786;
    B_97_321_0_idx_1 = rtb_B_97_30_0[1] * muDoubleScalarHypot(_rtB->B_97_307_0[1], _rtB->B_97_319_0[1]) * _rtP->P_783;
    B_5_0_0_idx_3 = (_rtP->P_784 * B_5_0_0_idx_3 - _rtP->P_785 * muDoubleScalarAtan2(_rtB->B_97_319_0[1], _rtB->B_97_307_0[1])) * _rtP->P_786;
    B_97_321_0_idx_2 = rtb_B_97_30_0[2] * muDoubleScalarHypot(_rtB->B_97_307_0[2], _rtB->B_97_319_0[2]) * _rtP->P_783;
    rtb_B_97_82_0 = (_rtP->P_784 * B_97_291_1_idx_2 - _rtP->P_785 * muDoubleScalarAtan2(_rtB->B_97_319_0[2], _rtB->B_97_307_0[2])) * _rtP->P_786;
    muDoubleScalarSinCos(rtb_B_97_66_0, &B_97_328_0_idx_0, &B_97_291_1_idx_2);
    muDoubleScalarSinCos(B_5_0_0_idx_3, &B_97_328_0_idx_1, &B_97_328_1_idx_1);
    muDoubleScalarSinCos(rtb_B_97_82_0, &B_97_328_0_idx_2, &B_97_328_1_idx_2);
    B_97_291_1_idx_2 *= B_97_321_0_idx_0;
    B_97_328_1_idx_1 *= B_97_321_0_idx_1;
    B_97_328_1_idx_2 *= B_97_321_0_idx_2;
    _rtB->B_97_330_0 = (B_97_291_1_idx_2 + B_97_328_1_idx_1) + B_97_328_1_idx_2;
}
_rtB->B_97_331_0 = _rtX->integ1_CSTATE;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T_PWORK.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T_PWORK.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_788;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        _rtB->B_97_332_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T_IWORK.CircularBufSize,
        &_rtDW->T_IWORK.Last,
        _rtDW->T_IWORK.Tail,
        _rtDW->T_IWORK.Head,
        _rtP->P_789,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
_rtB->B_97_333_0 = _rtB->B_97_331_0 - _rtB->B_97_332_0;
_rtB->B_97_334_0 = _rtX->Integ2_CSTATE;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_791;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        _rtB->B_97_335_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T1_IWORK.CircularBufSize,
        &_rtDW->T1_IWORK.Last,
        _rtDW->T1_IWORK.Tail,
        _rtDW->T1_IWORK.Head,
        _rtP->P_792,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
_rtB->B_97_336_0 = _rtB->B_97_334_0 - _rtB->B_97_335_0;
_rtB->B_97_337_0.re = _rtB->B_97_333_0;
_rtB->B_97_337_0.im = _rtB->B_97_336_0;
_rtB->B_97_338_0 = _rtX->integ1_CSTATE_n;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T_PWORK_l.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T_PWORK_l.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_794;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        _rtB->B_97_339_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T_IWORK_h.CircularBufSize,
        &_rtDW->T_IWORK_h.Last,
        _rtDW->T_IWORK_h.Tail,
        _rtDW->T_IWORK_h.Head,
        _rtP->P_795,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
_rtB->B_97_340_0 = _rtB->B_97_338_0 - _rtB->B_97_339_0;
_rtB->B_97_341_0 = _rtX->Integ2_CSTATE_p;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK_b.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK_b.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_797;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        _rtB->B_97_342_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T1_IWORK_f.CircularBufSize,
        &_rtDW->T1_IWORK_f.Last,
        _rtDW->T1_IWORK_f.Tail,
        _rtDW->T1_IWORK_f.Head,
        _rtP->P_798,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
_rtB->B_97_343_0 = _rtB->B_97_341_0 - _rtB->B_97_342_0;
_rtB->B_97_344_0.re = _rtB->B_97_340_0;
_rtB->B_97_344_0.im = _rtB->B_97_343_0;
_rtB->B_97_345_0 = _rtX->integ1_CSTATE_o;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T_PWORK_b.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T_PWORK_b.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_800;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        _rtB->B_97_346_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T_IWORK_n.CircularBufSize,
        &_rtDW->T_IWORK_n.Last,
        _rtDW->T_IWORK_n.Tail,
        _rtDW->T_IWORK_n.Head,
        _rtP->P_801,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
_rtB->B_97_347_0 = _rtB->B_97_345_0 - _rtB->B_97_346_0;
_rtB->B_97_348_0 = _rtX->Integ2_CSTATE_h;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK_j.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK_j.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_803;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        _rtB->B_97_349_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T1_IWORK_fz.CircularBufSize,
        &_rtDW->T1_IWORK_fz.Last,
        _rtDW->T1_IWORK_fz.Tail,
        _rtDW->T1_IWORK_fz.Head,
        _rtP->P_804,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
_rtB->B_97_350_0 = _rtB->B_97_348_0 - _rtB->B_97_349_0;
_rtB->B_97_351_0.re = _rtB->B_97_347_0;
_rtB->B_97_351_0.im = _rtB->B_97_350_0;
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    _rtB->B_97_353_0 = _rtB->B_97_352_0;
}


      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_NegSeqComputation(S, _rtB->B_97_353_0, _rtB->B_97_337_0, _rtB->B_97_344_0, _rtB->B_97_351_0, &_rtB->PosSeqComputation, &_rtDW->PosSeqComputation, &_rtP->PosSeqComputation);
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    _rtB->B_97_356_0 = _rtB->B_97_355_0;
}


      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_NegSeqComputation(S, _rtB->B_97_356_0, _rtB->B_97_337_0, _rtB->B_97_344_0, _rtB->B_97_351_0, &_rtB->NegSeqComputation, &_rtDW->NegSeqComputation, &_rtP->NegSeqComputation);
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    _rtB->B_97_359_0 = _rtB->B_97_358_0;
}


      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_ZeroSeqComputation(S, _rtB->B_97_359_0, _rtB->B_97_337_0, _rtB->B_97_344_0, _rtB->B_97_351_0, &_rtB->ZeroSeqComputation, &_rtDW->ZeroSeqComputation, &_rtP->ZeroSeqComputation);
_rtB->B_97_361_0[0] = muDoubleScalarHypot(_rtB->PosSeqComputation.B_39_2_0.re, _rtB->PosSeqComputation.B_39_2_0.im);
_rtB->B_97_361_0[1] = muDoubleScalarHypot(_rtB->NegSeqComputation.B_39_2_0.re, _rtB->NegSeqComputation.B_39_2_0.im);
_rtB->B_97_361_0[2] = muDoubleScalarHypot(_rtB->ZeroSeqComputation.B_41_1_0.re, _rtB->ZeroSeqComputation.B_41_1_0.im);
rtb_B_97_30_0[0] = muDoubleScalarAtan2(_rtB->PosSeqComputation.B_39_2_0.im, _rtB->PosSeqComputation.B_39_2_0.re);
rtb_B_97_30_0[1] = muDoubleScalarAtan2(_rtB->NegSeqComputation.B_39_2_0.im, _rtB->NegSeqComputation.B_39_2_0.re);
rtb_B_97_30_0[2] = muDoubleScalarAtan2(_rtB->ZeroSeqComputation.B_41_1_0.im, _rtB->ZeroSeqComputation.B_41_1_0.re);
_rtB->B_97_362_0 = _rtP->P_808 * _rtB->B_97_361_0[0];
_rtB->B_97_363_0 = 0.0;
_rtB->B_97_363_0 += _rtP->P_810 * _rtX->TransferFcn1_CSTATE;
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    B_97_321_0_idx_0 *= B_97_328_0_idx_0;
    B_97_321_0_idx_1 *= B_97_328_0_idx_1;
    B_97_321_0_idx_2 *= B_97_328_0_idx_2;
    _rtB->B_97_365_0 = (B_97_321_0_idx_0 + B_97_321_0_idx_1) + B_97_321_0_idx_2;
}
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 366, SS_CALL_MDL_OUTPUTS);
   
  
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    for (i = 0; i < 6; i++) {
        if (_rtB->B_97_116_1[i + 13] >= _rtP->P_813) {
            B_5_0_0_idx_3 = _rtB->B_97_116_0[i + 13] * _rtP->P_812;
        } else {
            B_5_0_0_idx_3 = _rtP->P_811;
        }
        if (B_5_0_0_idx_3 > _rtP->P_814) {
            _rtB->B_97_370_0[i] = _rtP->P_814;
        } else if (B_5_0_0_idx_3 < _rtP->P_815) {
            _rtB->B_97_370_0[i] = _rtP->P_815;
        } else {
            _rtB->B_97_370_0[i] = B_5_0_0_idx_3;
        }
    }
    for (i = 0; i < 6; i++) {
        if (_rtDW->UnitDelay_DSTATE_ew[i] != 0.0) {
            B_97_374_0[i] = _rtB->B_97_372_0[i];
        } else {
            B_97_374_0[i] = _rtB->B_97_373_0[i];
        }
    }
}
i = ssIsSampleHit(S, 3, 0);
if (i != 0) {
    _rtB->B_97_375_0[0] = _rtDW->UnitDelay2_DSTATE[0];
    _rtB->B_97_375_0[1] = _rtDW->UnitDelay2_DSTATE[1];
    _rtB->B_97_375_0[2] = _rtDW->UnitDelay2_DSTATE[2];
}
_rtB->B_97_382_0 = look1_binlxpw(muDoubleScalarRem(ssGetT(S) - _rtB->B_97_377_0, _rtB->B_97_379_0), _rtP->P_822, _rtP->P_821, 3U);
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    i = (_rtB->B_97_375_0[0] >= _rtB->B_97_382_0);
    rtb_B_97_384_0 = (i != 0);
    rtb_B_97_76_0 = !rtb_B_97_384_0;
    B_97_291_1_idx_2 = i;
    rtb_B_97_111_0 = rtb_B_97_384_0;
    i = (_rtB->B_97_375_0[1] >= _rtB->B_97_382_0);
    rtb_B_97_384_0 = (i != 0);
    B_97_385_0_idx_1 = !rtb_B_97_384_0;
    B_97_328_1_idx_1 = i;
    B_97_384_0_idx_1 = rtb_B_97_384_0;
    i = (_rtB->B_97_375_0[2] >= _rtB->B_97_382_0);
    rtb_B_97_384_0 = (i != 0);
    B_97_328_1_idx_2 = i;
    _rtB->B_97_387_0[0] = rtb_B_97_111_0;
    _rtB->B_97_387_0[1] = rtb_B_97_76_0;
    _rtB->B_97_387_0[2] = B_97_384_0_idx_1;
    _rtB->B_97_387_0[3] = B_97_385_0_idx_1;
    _rtB->B_97_387_0[4] = rtb_B_97_384_0;
    _rtB->B_97_387_0[5] = !rtb_B_97_384_0;
    tmpForInput[0] = B_97_212_0_idx_0;
    tmpForInput[1] = B_97_212_0_idx_1;
    tmpForInput[2] = B_97_212_0_idx_2;
    tmpForInput[3] = rtb_B_97_215_0;
    tmpForInput[4] = rtb_B_97_216_0;
    tmpForInput[5] = rtb_B_97_221_0;
    tmpForInput[6] = rtb_B_97_224_0;
    tmpForInput[7] = rtb_B_97_225_0;
    tmpForInput[8] = rtb_B_97_226_0;
    rtb_B_97_226_0 = -0.0;
    for (i = 0; i < 9; i++) {
        rtb_B_97_226_0 += tmpForInput[i];
    }
    B_97_395_0 = _rtP->P_823 * rtb_B_97_226_0;
    tmpForInput[0] = B_97_173_0_idx_0;
    tmpForInput[1] = B_97_173_0_idx_1;
    tmpForInput[2] = B_97_173_0_idx_2;
    tmpForInput[3] = rtb_B_97_176_0;
    tmpForInput[4] = rtb_B_97_177_0;
    tmpForInput[5] = rtb_B_97_182_0;
    tmpForInput[6] = rtb_B_97_185_0;
    tmpForInput[7] = rtb_B_97_186_0;
    tmpForInput[8] = rtb_B_97_187_0;
    rtb_B_97_396_0 = -0.0;
    for (i = 0; i < 9; i++) {
        rtb_B_97_396_0 += tmpForInput[i];
    }
    B_97_397_0 = _rtP->P_824 * rtb_B_97_396_0;
}
if (ssGetT(S) >= _rtP->P_829) {
    _rtB->B_97_405_0 = _rtB->B_97_400_0;
} else {
    if (ssGetT(S) >= _rtP->P_500) {
        _rtB->B_88_0_0 = _rtB->B_97_402_0;
    } else {
        _rtB->B_88_0_0 = _rtB->B_97_400_0;
    }
    _rtB->B_97_405_0 = _rtB->B_88_0_0;
}
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->TransportDelay1_PWORK.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->TransportDelay1_PWORK.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_830;    
    
       
    
    
  
    
      
                     


          
        

      
      
          if (_rtP->P_830 == 0.0)
          _rtB->B_97_406_0 = _rtB->B_97_405_0; 
          else
        
        _rtB->B_97_406_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->TransportDelay1_IWORK.CircularBufSize,
        &_rtDW->TransportDelay1_IWORK.Last,
        _rtDW->TransportDelay1_IWORK.Tail,
        _rtDW->TransportDelay1_IWORK.Head,
        _rtP->P_831,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
_rtB->B_97_407_0 = _rtB->B_97_406_0;
i = ssIsSampleHit(S, 1, 0);
if ((i != 0) && (ssIsMajorTimeStep(S) != 0)) {
    if (_rtB->B_97_407_0 > 0.0) {
        if (!_rtDW->GridSupportingasCurrentSourceGridFeeding_MODE) {
            if (ssGetTaskTime(S, 1) != ssGetTStart(S)) {
                ssSetBlockStateForSolverChangedAtMajorStep(S);
            }
    


      
      (void) memset(&(((XDis_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetContStateDisabled(S))->TransferFcn1_CSTATE_p), 0,
17*sizeof(boolean_T));

            _rtDW->Integ4_SYSTEM_ENABLE_fj = 1U;
            _rtDW->GridSupportingasCurrentSourceGridFeeding_MODE = true;
        }
    } else {
        if (_rtDW->GridSupportingasCurrentSourceGridFeeding_MODE) {
            ssSetBlockStateForSolverChangedAtMajorStep(S);
    


      
      (void) memset(&(((XDis_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetContStateDisabled(S))->TransferFcn1_CSTATE_p), 1,
17*sizeof(boolean_T));

            if (_rtDW->AutomaticGainControl_MODE_c) {
                _rtDW->AutomaticGainControl_MODE_c = false;
            }
            _rtDW->GridSupportingasCurrentSourceGridFeeding_MODE = false;
        }
    }
}
if (_rtDW->GridSupportingasCurrentSourceGridFeeding_MODE) {
    _rtB->B_27_0_0 = 0.0;
    _rtB->B_27_0_0 += _rtP->P_59 * _rtX->TransferFcn1_CSTATE_p;
    _rtB->B_27_2_0 = _rtB->B_27_0_0 - _rtB->B_27_1_0;
    if (ssGetT(S) >= _rtP->P_62) {
        _rtB->B_27_5_0 = _rtB->B_27_2_0;
    } else {
        _rtB->B_27_5_0 = _rtB->B_27_4_0;
    }
    _rtB->B_27_6_0 = _rtX->Integrator_CSTATE_e;
    _rtB->B_27_7_0 = _rtP->P_64 * _rtB->B_27_5_0;
    _rtB->B_27_8_0 = _rtX->Filter_CSTATE_pr;
    _rtB->B_27_9_0 = _rtB->B_27_7_0 - _rtB->B_27_8_0;
    _rtB->B_27_10_0 = _rtP->P_66 * _rtB->B_27_9_0;
    _rtB->B_27_11_0 = 0.0;
    _rtB->B_27_11_0 += _rtP->P_68 * _rtX->TransferFcn2_CSTATE_f;
    _rtB->B_27_13_0 = _rtX->Integrator_CSTATE_o;
    _rtB->B_27_14_0 = _rtP->P_71 * _rtB->B_27_11_0;
    _rtB->B_27_15_0 = _rtX->Filter_CSTATE_d;
    _rtB->B_27_16_0 = _rtB->B_27_14_0 - _rtB->B_27_15_0;
    _rtB->B_27_17_0 = _rtP->P_73 * _rtB->B_27_16_0;
    _rtB->B_27_19_0 = 0.0;
    _rtB->B_27_19_0 += _rtP->P_75 * _rtX->TransferFcn4_CSTATE_p;
    _rtB->B_27_20_0 = _rtX->Integrator_CSTATE_p;
    _rtB->B_27_21_0 = _rtP->P_77 * _rtB->B_27_19_0;
    _rtB->B_27_22_0 = _rtX->Filter_CSTATE_j;
    _rtB->B_27_23_0 = _rtB->B_27_21_0 - _rtB->B_27_22_0;
    _rtB->B_27_24_0 = _rtP->P_79 * _rtB->B_27_23_0;
    if (ssGetT(S) >= _rtP->P_80) {
        rtb_B_97_176_0 = ((_rtP->P_23 * _rtB->B_27_5_0 + _rtB->B_27_6_0) + _rtB->B_27_10_0) + ((_rtP->P_69 * _rtB->B_27_11_0 + _rtB->B_27_13_0) + _rtB->B_27_17_0);
    } else {
        rtb_B_97_176_0 = (_rtP->P_69 * _rtB->B_27_11_0 + _rtB->B_27_13_0) + _rtB->B_27_17_0;
    }
    _rtB->B_27_31_0 = _rtP->P_81 * _rtB->B_27_11_0;
    _rtB->B_27_32_0 = _rtP->P_82 * _rtB->B_27_5_0;
    _rtB->B_27_33_0 = _rtP->P_83 * _rtB->B_27_19_0;
    if (ssGetT(S) >= _rtP->P_85) {
        _rtB->B_27_35_0 = _rtB->B_27_2_0;
    } else {
        _rtB->B_27_35_0 = _rtB->B_27_34_0;
    }
    if (ssGetT(S) >= _rtP->P_87) {
        _rtB->B_27_37_0 = _rtB->B_27_2_0;
    } else {
        _rtB->B_27_37_0 = _rtB->B_27_36_0;
    }
    if (ssGetT(S) >= _rtP->P_88) {
        rtb_B_97_176_0 += (_rtP->P_22 * _rtB->B_27_19_0 + _rtB->B_27_20_0) + _rtB->B_27_24_0;
    }
    _rtB->B_27_40_0 = 0.0;
    _rtB->B_27_40_0 += _rtP->P_90 * _rtX->TransferFcn4_CSTATE_d;
    _rtB->B_27_41_0 = _rtX->Integrator_CSTATE_pv;
    _rtB->B_27_42_0 = _rtP->P_92 * _rtB->B_27_40_0;
    _rtB->B_27_43_0 = _rtX->Filter_CSTATE_l;
    _rtB->B_27_44_0 = _rtB->B_27_42_0 - _rtB->B_27_43_0;
    _rtB->B_27_45_0 = _rtP->P_94 * _rtB->B_27_44_0;
    _rtB->B_27_47_0 = _rtB->B_97_362_0 - _rtB->B_27_46_0;
    if (ssGetT(S) >= _rtP->P_97) {
        _rtB->B_27_50_0 = _rtB->B_27_47_0;
    } else {
        _rtB->B_27_50_0 = _rtB->B_27_49_0;
    }
    _rtB->B_27_51_0 = _rtX->Integrator_CSTATE_ou;
    _rtB->B_27_52_0 = _rtP->P_99 * _rtB->B_27_50_0;
    _rtB->B_27_53_0 = _rtX->Filter_CSTATE_f;
    _rtB->B_27_54_0 = _rtB->B_27_52_0 - _rtB->B_27_53_0;
    _rtB->B_27_55_0 = _rtP->P_101 * _rtB->B_27_54_0;
    _rtB->B_27_56_0 = 0.0;
    _rtB->B_27_56_0 += _rtP->P_103 * _rtX->TransferFcn2_CSTATE_b;
    _rtB->B_27_58_0 = _rtX->Integrator_CSTATE_cm;
    _rtB->B_27_59_0 = _rtP->P_106 * _rtB->B_27_56_0;
    _rtB->B_27_60_0 = _rtX->Filter_CSTATE_n;
    _rtB->B_27_61_0 = _rtB->B_27_59_0 - _rtB->B_27_60_0;
    _rtB->B_27_62_0 = _rtP->P_108 * _rtB->B_27_61_0;
    if (ssGetT(S) >= _rtP->P_109) {
        rtb_B_97_82_0 = ((_rtP->P_57 * _rtB->B_27_50_0 + _rtB->B_27_51_0) + _rtB->B_27_55_0) + ((_rtP->P_104 * _rtB->B_27_56_0 + _rtB->B_27_58_0) + _rtB->B_27_62_0);
    } else {
        rtb_B_97_82_0 = (_rtP->P_104 * _rtB->B_27_56_0 + _rtB->B_27_58_0) + _rtB->B_27_62_0;
    }
    if (ssGetT(S) >= _rtP->P_110) {
        rtb_B_97_82_0 += (_rtP->P_56 * _rtB->B_27_40_0 + _rtB->B_27_41_0) + _rtB->B_27_45_0;
    }
    _rtB->B_27_72_0 = _rtP->P_112 * rtb_B_97_176_0 + _rtB->B_97_398_0;
    if (ssIsMajorTimeStep(S) != 0) {
        _rtDW->Saturation_MODE_i = _rtB->B_27_72_0 >= _rtP->P_113 ? 1 : _rtB->B_27_72_0 > _rtP->P_114 ? 0 : -1;
    }
    rtb_B_97_176_0 = _rtDW->Saturation_MODE_i == 1 ? _rtP->P_113 : _rtDW->Saturation_MODE_i == -1 ? _rtP->P_114 : _rtB->B_27_72_0;
    _rtB->B_27_74_0 = _rtP->P_111 * rtb_B_97_82_0 + _rtB->B_97_399_0;
    if (ssIsMajorTimeStep(S) != 0) {
        _rtDW->Saturation1_MODE_d = _rtB->B_27_74_0 >= _rtP->P_115 ? 1 : _rtB->B_27_74_0 > _rtP->P_116 ? 0 : -1;
    }
    rtb_B_97_82_0 = _rtDW->Saturation1_MODE_d == 1 ? _rtP->P_115 : _rtDW->Saturation1_MODE_d == -1 ? _rtP->P_116 : _rtB->B_27_74_0;
    i = ssIsSampleHit(S, 2, 0);
    if (i != 0) {
        _rtB->B_27_78_0 = muDoubleScalarMod(_rtDW->DiscreteTimeIntegrator_DSTATE_i, _rtB->B_27_77_0);
        B_27_79_0_idx_0 = _rtP->P_120 * B_97_173_0_idx_0;
        B_27_79_0_idx_1 = _rtP->P_120 * B_97_173_0_idx_1;
        B_27_79_0_idx_2 = _rtP->P_120 * B_97_173_0_idx_2;
        rtb_B_27_80_0 = muDoubleScalarSin(_rtB->B_27_78_0);
        rtb_B_27_81_0 = muDoubleScalarCos(_rtB->B_27_78_0);
        rtb_B_27_86_0 = (0.0 - rtb_B_27_80_0 * _rtB->B_27_82_0) - rtb_B_27_81_0 * _rtB->B_27_84_0;
        rtb_B_27_89_0 = rtb_B_27_80_0 * _rtB->B_27_84_0 - rtb_B_27_81_0 * _rtB->B_27_82_0;
        rtb_B_27_90_0 = (0.0 - rtb_B_27_86_0) - rtb_B_27_80_0;
        rtb_B_27_91_0 = (0.0 - rtb_B_27_89_0) - rtb_B_27_81_0;
        _rtB->B_27_92_0 = ((B_27_79_0_idx_0 * rtb_B_27_80_0 + B_27_79_0_idx_1 * rtb_B_27_86_0) + B_27_79_0_idx_2 * rtb_B_27_90_0) * 0.66666666666666663;
        _rtB->B_27_93_0 = ((B_27_79_0_idx_0 * rtb_B_27_81_0 + B_27_79_0_idx_1 * rtb_B_27_89_0) + B_27_79_0_idx_2 * rtb_B_27_91_0) * 0.66666666666666663;
        B_27_95_0_idx_0 = _rtP->P_123 * B_97_212_0_idx_0 * _rtP->P_124;
        B_27_95_0_idx_1 = _rtP->P_123 * B_97_212_0_idx_1 * _rtP->P_124;
        B_27_95_0_idx_2 = _rtP->P_123 * B_97_212_0_idx_2 * _rtP->P_124;
        rtb_B_27_100_0 = (0.0 - rtb_B_27_80_0 * _rtB->B_27_96_0) - rtb_B_27_81_0 * _rtB->B_27_98_0;
        rtb_B_27_103_0 = rtb_B_27_80_0 * _rtB->B_27_98_0 - rtb_B_27_81_0 * _rtB->B_27_96_0;
        rtb_B_27_104_0 = (0.0 - rtb_B_27_100_0) - rtb_B_27_80_0;
        rtb_B_27_105_0 = (0.0 - rtb_B_27_103_0) - rtb_B_27_81_0;
        _rtB->B_27_106_0 = ((B_27_95_0_idx_0 * rtb_B_27_80_0 + B_27_95_0_idx_1 * rtb_B_27_100_0) + B_27_95_0_idx_2 * rtb_B_27_104_0) * 0.66666666666666663;
        _rtB->B_27_107_0 = ((B_27_95_0_idx_0 * rtb_B_27_81_0 + B_27_95_0_idx_1 * rtb_B_27_103_0) + B_27_95_0_idx_2 * rtb_B_27_105_0) * 0.66666666666666663;
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 27, 108, SS_CALL_MDL_OUTPUTS);
   
  
    }
    _rtB->B_27_110_0 = (_rtB->B_97_259_0[0] - rtb_B_97_176_0) * _rtP->P_127;
    i = ssIsSampleHit(S, 2, 0);
    if (i != 0) {
        _rtB->B_27_112_0 = _rtP->P_129 * _rtB->B_27_110_0;
        _rtB->B_27_113_0 = _rtP->P_130 * _rtB->B_27_112_0 + _rtDW->Integrator_DSTATE_b;
        B_5_0_0_idx_3 = _rtP->P_128 * _rtB->B_27_110_0 + _rtB->B_27_113_0;
        if (B_5_0_0_idx_3 > _rtP->P_132) {
            _rtB->B_27_115_0 = _rtP->P_132;
        } else if (B_5_0_0_idx_3 < _rtP->P_133) {
            _rtB->B_27_115_0 = _rtP->P_133;
        } else {
            _rtB->B_27_115_0 = B_5_0_0_idx_3;
        }
    }
    _rtB->B_27_117_0 = (_rtB->B_97_259_0[1] - rtb_B_97_82_0) * _rtP->P_134;
    i = ssIsSampleHit(S, 2, 0);
    if (i != 0) {
        _rtB->B_27_119_0 = _rtP->P_136 * _rtB->B_27_117_0;
        _rtB->B_27_120_0 = _rtP->P_137 * _rtB->B_27_119_0 + _rtDW->Integrator_DSTATE_e;
        B_5_0_0_idx_3 = _rtP->P_135 * _rtB->B_27_117_0 + _rtB->B_27_120_0;
        if (B_5_0_0_idx_3 > _rtP->P_139) {
            _rtB->B_27_122_0 = _rtP->P_139;
        } else if (B_5_0_0_idx_3 < _rtP->P_140) {
            _rtB->B_27_122_0 = _rtP->P_140;
        } else {
            _rtB->B_27_122_0 = B_5_0_0_idx_3;
        }
        _rtB->B_27_123_0[0] = _rtB->B_27_115_0 - _rtB->B_27_106_0;
        _rtB->B_27_123_0[1] = _rtB->B_27_122_0 - _rtB->B_27_107_0;
        B_5_0_0_idx_3 = _rtP->P_141 * _rtB->B_27_123_0[0] + _rtDW->Integrator_DSTATE_b2[0];
        if (B_5_0_0_idx_3 > _rtP->P_144) {
            _rtB->B_27_127_0[0] = _rtP->P_144;
        } else if (B_5_0_0_idx_3 < _rtP->P_145) {
            _rtB->B_27_127_0[0] = _rtP->P_145;
        } else {
            _rtB->B_27_127_0[0] = B_5_0_0_idx_3;
        }
        B_5_0_0_idx_3 = _rtP->P_141 * _rtB->B_27_123_0[1] + _rtDW->Integrator_DSTATE_b2[1];
        if (B_5_0_0_idx_3 > _rtP->P_144) {
            _rtB->B_27_127_0[1] = _rtP->P_144;
        } else if (B_5_0_0_idx_3 < _rtP->P_145) {
            _rtB->B_27_127_0[1] = _rtP->P_145;
        } else {
            _rtB->B_27_127_0[1] = B_5_0_0_idx_3;
        }
        _rtB->B_27_130_0 = (_rtP->P_146 * _rtB->B_27_115_0 + _rtB->B_27_92_0) - _rtP->P_147 * _rtB->B_27_122_0;
        _rtB->B_27_133_0 = (_rtP->P_148 * _rtB->B_27_122_0 + _rtB->B_27_93_0) + _rtP->P_149 * _rtB->B_27_115_0;
        _rtB->B_27_135_0[0] = (_rtB->B_27_127_0[0] + _rtB->B_27_130_0) * _rtP->P_150;
        _rtB->B_27_135_0[1] = (_rtB->B_27_127_0[1] + _rtB->B_27_133_0) * _rtP->P_150;
        if (_rtB->B_27_135_0[0] > _rtP->P_151) {
            _rtB->B_27_136_0[0] = _rtP->P_151;
        } else if (_rtB->B_27_135_0[0] < _rtP->P_152) {
            _rtB->B_27_136_0[0] = _rtP->P_152;
        } else {
            _rtB->B_27_136_0[0] = _rtB->B_27_135_0[0];
        }
        if (_rtB->B_27_135_0[1] > _rtP->P_151) {
            _rtB->B_27_136_0[1] = _rtP->P_151;
        } else if (_rtB->B_27_135_0[1] < _rtP->P_152) {
            _rtB->B_27_136_0[1] = _rtP->P_152;
        } else {
            _rtB->B_27_136_0[1] = _rtB->B_27_135_0[1];
        }
        B_5_0_0_idx_3 = _rtB->B_27_136_0[0] / _rtB->B_27_140_0;
        rtb_B_97_176_0 = _rtB->B_27_136_0[1] / _rtB->B_27_140_0;
        rtb_B_97_82_0 = muDoubleScalarAtan2(rtb_B_97_176_0, B_5_0_0_idx_3);
        rtb_B_97_176_0 = muDoubleScalarHypot(B_5_0_0_idx_3, rtb_B_97_176_0);
        if (rtb_B_97_176_0 > _rtP->P_156) {
            rtb_B_97_176_0 = _rtP->P_156;
        } else {
            if (rtb_B_97_176_0 < _rtP->P_157) {
                rtb_B_97_176_0 = _rtP->P_157;
            }
        }
        _rtB->B_27_151_0[0] = muDoubleScalarSin((((_rtB->B_27_78_0 + _rtB->B_27_145_0[0]) + _rtB->B_27_146_0) + _rtB->B_27_147_0) + rtb_B_97_82_0) * rtb_B_97_176_0;
        _rtB->B_27_151_0[1] = muDoubleScalarSin((((_rtB->B_27_78_0 + _rtB->B_27_145_0[1]) + _rtB->B_27_146_0) + _rtB->B_27_147_0) + rtb_B_97_82_0) * rtb_B_97_176_0;
        _rtB->B_27_151_0[2] = muDoubleScalarSin((((_rtB->B_27_78_0 + _rtB->B_27_145_0[2]) + _rtB->B_27_146_0) + _rtB->B_27_147_0) + rtb_B_97_82_0) * rtb_B_97_176_0;
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 27, 152, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 27, 153, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 27, 154, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 27, 155, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 27, 156, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 27, 157, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 27, 158, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 27, 159, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 27, 160, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 27, 161, SS_CALL_MDL_OUTPUTS);
   
  
        _rtB->B_27_162_0[0] = _rtP->P_161 * _rtB->B_27_123_0[0];
        _rtB->B_27_162_0[1] = _rtP->P_161 * _rtB->B_27_123_0[1];
    }
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 27, 163, SS_CALL_MDL_OUTPUTS);
   
  
    i = ssIsSampleHit(S, 2, 0);
    if (i != 0) {
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 27, 164, SS_CALL_MDL_OUTPUTS);
   
  
        rtb_B_97_176_0 = (0.0 - rtb_B_27_80_0 * _rtB->B_27_165_0) - rtb_B_27_81_0 * _rtB->B_27_167_0;
        rtb_B_97_226_0 = rtb_B_27_80_0 * _rtB->B_27_167_0 - rtb_B_27_81_0 * _rtB->B_27_165_0;
        rtb_B_97_396_0 = (0.0 - rtb_B_97_176_0) - rtb_B_27_80_0;
        rtb_B_97_177_0 = (0.0 - rtb_B_97_226_0) - rtb_B_27_81_0;
        _rtB->B_27_175_0 = ((B_27_79_0_idx_0 * rtb_B_27_81_0 + B_27_79_0_idx_1 * rtb_B_97_226_0) + B_27_79_0_idx_2 * rtb_B_97_177_0) * 0.66666666666666663;
        if (_rtDW->Integ4_SYSTEM_ENABLE_fj != 0) {
            _rtB->B_27_176_0 = _rtDW->Integ4_DSTATE_p3;
        } else {
            _rtB->B_27_176_0 = _rtP->P_164 * _rtB->B_27_175_0 + _rtDW->Integ4_DSTATE_p3;
        }
        rtb_B_97_182_0 = _rtDW->UnitDelay_DSTATE_oo;
        if (_rtDW->UnitDelay_DSTATE_oo > _rtP->P_167) {
            B_5_0_0_idx_3 = _rtP->P_167;
        } else if (_rtDW->UnitDelay_DSTATE_oo < _rtP->P_168) {
            B_5_0_0_idx_3 = _rtP->P_168;
        } else {
            B_5_0_0_idx_3 = _rtDW->UnitDelay_DSTATE_oo;
        }
        rtb_B_97_185_0 = 1.0 / B_5_0_0_idx_3 / 1.0e-5;
        rtb_B_97_186_0 = muDoubleScalarCeil(rtb_B_97_185_0);
        _rtB->B_27_181_0 = _rtP->P_169 * rtb_B_97_186_0;
                
          /* Level2 S-Function Block: '<S484>/B_27_17' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 27, 182, SS_CALL_MDL_OUTPUTS);

   
  
        rtb_B_97_82_0 = ssGetTaskTime(S, 2);
        if (rtb_B_97_82_0 < _rtP->P_171) {
            rtb_B_97_82_0 = _rtP->P_172;
        } else {
            rtb_B_97_82_0 = _rtP->P_173;
        }
        if (rtb_B_97_82_0 >= _rtP->P_175) {
            rtb_B_97_82_0 = rtb_B_97_185_0 - rtb_B_97_186_0;
            rtb_B_97_185_0 = ((_rtB->B_27_175_0 - _rtDW->UnitDelay_DSTATE_k2) * rtb_B_97_82_0 * _rtP->P_55 + _rtB->B_27_175_0) * (rtb_B_97_82_0 / rtb_B_97_185_0) + (_rtB->B_27_176_0 - _rtB->B_27_182_0) * rtb_B_97_182_0;
        } else {
            rtb_B_97_185_0 = _rtP->P_174;
        }
        _rtB->B_27_189_0 = _rtB->B_27_188_0;
        i = ssIsSampleHit(S, 2, 0);
        if ((i != 0) && (ssIsMajorTimeStep(S) != 0)) {
            if (_rtB->B_27_189_0 > 0.0) {
                if (!_rtDW->AutomaticGainControl_MODE_c) {
                    if (ssGetTaskTime(S, 2) != ssGetTStart(S)) {
                        ssSetBlockStateForSolverChangedAtMajorStep(S);
                    }
  
                    _rtDW->Integ4_SYSTEM_ENABLE_gq = 1U;
                    _rtDW->Integ4_SYSTEM_ENABLE_li = 1U;
                    _rtDW->AutomaticGainControl_MODE_c = true;
                }
            } else {
                if (_rtDW->AutomaticGainControl_MODE_c) {
                    ssSetBlockStateForSolverChangedAtMajorStep(S);
  
                    _rtDW->AutomaticGainControl_MODE_c = false;
                }
            }
        }
        if (_rtDW->AutomaticGainControl_MODE_c) {
            i = ssIsSampleHit(S, 2, 0);
            if (i != 0) {
                rtb_B_97_186_0 = (0.0 - rtb_B_27_80_0 * _rtB->B_23_0_0) - rtb_B_27_81_0 * _rtB->B_23_2_0;
                rtb_B_97_187_0 = rtb_B_27_80_0 * _rtB->B_23_2_0 - rtb_B_27_81_0 * _rtB->B_23_0_0;
                rtb_B_97_215_0 = (0.0 - rtb_B_97_186_0) - rtb_B_27_80_0;
                rtb_B_97_216_0 = (0.0 - rtb_B_97_187_0) - rtb_B_27_81_0;
                _rtB->B_23_10_0 = ((B_27_79_0_idx_0 * rtb_B_27_80_0 + B_27_79_0_idx_1 * rtb_B_97_186_0) + B_27_79_0_idx_2 * rtb_B_97_215_0) * 0.66666666666666663;
                if (_rtDW->Integ4_SYSTEM_ENABLE_gq != 0) {
                    _rtB->B_23_11_0 = _rtDW->Integ4_DSTATE_jl;
                } else {
                    _rtB->B_23_11_0 = _rtP->P_29 * _rtB->B_23_10_0 + _rtDW->Integ4_DSTATE_jl;
                }
                if (rtb_B_97_182_0 > _rtP->P_31) {
                    rtb_B_97_82_0 = _rtP->P_31;
                } else if (rtb_B_97_182_0 < _rtP->P_32) {
                    rtb_B_97_82_0 = _rtP->P_32;
                } else {
                    rtb_B_97_82_0 = rtb_B_97_182_0;
                }
                rtb_B_97_221_0 = 1.0 / rtb_B_97_82_0 / 1.0e-5;
                rtb_B_97_224_0 = muDoubleScalarCeil(rtb_B_97_221_0);
                _rtB->B_23_15_0 = _rtP->P_33 * rtb_B_97_224_0;
                
          /* Level2 S-Function Block: '<S479>/B_23_0' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 23, 16, SS_CALL_MDL_OUTPUTS);

   
  
                rtb_B_97_82_0 = ssGetTaskTime(S, 2);
                if (rtb_B_97_82_0 < _rtP->P_35) {
                    rtb_B_97_82_0 = _rtP->P_36;
                } else {
                    rtb_B_97_82_0 = _rtP->P_37;
                }
                if (rtb_B_97_82_0 >= _rtP->P_39) {
                    rtb_B_97_82_0 = rtb_B_97_221_0 - rtb_B_97_224_0;
                    rtb_B_97_221_0 = ((_rtB->B_23_10_0 - _rtDW->UnitDelay_DSTATE_cu) * rtb_B_97_82_0 * _rtP->P_24 + _rtB->B_23_10_0) * (rtb_B_97_82_0 / rtb_B_97_221_0) + (_rtB->B_23_11_0 - _rtB->B_23_16_0) * rtb_B_97_182_0;
                } else {
                    rtb_B_97_221_0 = _rtP->P_38;
                }
                _rtB->B_23_22_0 = ((B_27_79_0_idx_0 * rtb_B_27_81_0 + B_27_79_0_idx_1 * rtb_B_97_187_0) + B_27_79_0_idx_2 * rtb_B_97_216_0) * 0.66666666666666663;
                if (_rtDW->Integ4_SYSTEM_ENABLE_li != 0) {
                    _rtB->B_23_23_0 = _rtDW->Integ4_DSTATE_it;
                } else {
                    _rtB->B_23_23_0 = _rtP->P_40 * _rtB->B_23_22_0 + _rtDW->Integ4_DSTATE_it;
                }
                if (rtb_B_97_182_0 > _rtP->P_42) {
                    rtb_B_97_82_0 = _rtP->P_42;
                } else if (rtb_B_97_182_0 < _rtP->P_43) {
                    rtb_B_97_82_0 = _rtP->P_43;
                } else {
                    rtb_B_97_82_0 = rtb_B_97_182_0;
                }
                rtb_B_97_224_0 = 1.0 / rtb_B_97_82_0 / 1.0e-5;
                rtb_B_97_225_0 = muDoubleScalarCeil(rtb_B_97_224_0);
                _rtB->B_23_27_0 = _rtP->P_44 * rtb_B_97_225_0;
                
          /* Level2 S-Function Block: '<S481>/B_23_2' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 23, 28, SS_CALL_MDL_OUTPUTS);

   
  
                rtb_B_97_82_0 = ssGetTaskTime(S, 2);
                if (rtb_B_97_82_0 < _rtP->P_46) {
                    rtb_B_97_82_0 = _rtP->P_47;
                } else {
                    rtb_B_97_82_0 = _rtP->P_48;
                }
                if (rtb_B_97_82_0 >= _rtP->P_50) {
                    rtb_B_97_82_0 = rtb_B_97_224_0 - rtb_B_97_225_0;
                    rtb_B_97_82_0 = ((_rtB->B_23_22_0 - _rtDW->UnitDelay_DSTATE_eb) * rtb_B_97_82_0 * _rtP->P_25 + _rtB->B_23_22_0) * (rtb_B_97_82_0 / rtb_B_97_224_0) + (_rtB->B_23_23_0 - _rtB->B_23_28_0) * rtb_B_97_182_0;
                } else {
                    rtb_B_97_82_0 = _rtP->P_49;
                }
                B_23_36_0 = _rtP->P_51 * muDoubleScalarAtan2(rtb_B_97_82_0, rtb_B_97_221_0);
                tmpForInput[0] = B_27_79_0_idx_0;
                tmpForInput[1] = B_27_79_0_idx_1;
                tmpForInput[2] = B_27_79_0_idx_2;
                tmpForInput[3] = rtb_B_27_80_0;
                tmpForInput[4] = rtb_B_27_81_0;
                tmpForInput[5] = rtb_B_97_186_0;
                tmpForInput[6] = rtb_B_97_187_0;
                tmpForInput[7] = rtb_B_97_215_0;
                tmpForInput[8] = rtb_B_97_216_0;
                B_5_0_0_idx_3 = -0.0;
                for (i = 0; i < 9; i++) {
                    B_5_0_0_idx_3 += tmpForInput[i];
                }
                B_23_38_0 = _rtP->P_52 * B_5_0_0_idx_3;
                B_5_0_0_idx_3 = muDoubleScalarHypot(rtb_B_97_221_0, rtb_B_97_82_0);
                if (B_5_0_0_idx_3 > _rtP->P_53) {
                    B_5_0_0_idx_3 = _rtP->P_53;
                } else {
                    if (B_5_0_0_idx_3 < _rtP->P_54) {
                        B_5_0_0_idx_3 = _rtP->P_54;
                    }
                }
                _rtB->B_23_40_0 = 1.0 / B_5_0_0_idx_3;
            }
            if (ssIsMajorTimeStep(S) != 0) {
                srUpdateBC(_rtDW->AutomaticGainControl_SubsysRanBC_nm);
            }
        }
        rtb_B_97_82_0 = rtb_B_97_185_0 * _rtB->B_23_40_0;
        _rtB->B_27_196_0 = _rtP->P_182 * rtb_B_97_82_0 * _rtP->P_183;
        B_27_197_0 = _rtDW->UD_DSTATE_i;
        B_27_198_0 = _rtB->B_27_196_0 - B_27_197_0;
        B_5_0_0_idx_3 = (_rtP->P_177 * rtb_B_97_82_0 + _rtDW->DiscreteTimeIntegrator_DSTATE_oj) + B_27_198_0;
        if (B_5_0_0_idx_3 > _rtP->P_185) {
            _rtB->B_27_200_0 = _rtP->P_185;
        } else if (B_5_0_0_idx_3 < _rtP->P_186) {
            _rtB->B_27_200_0 = _rtP->P_186;
        } else {
            _rtB->B_27_200_0 = B_5_0_0_idx_3;
        }
        B_5_0_0_idx_3 = _rtP->P_187 * _rtB->B_27_200_0 - _rtDW->UnitDelay_DSTATE_b4;
        if (B_5_0_0_idx_3 > _rtP->P_189) {
            B_5_0_0_idx_3 = _rtP->P_189;
        } else {
            if (B_5_0_0_idx_3 < _rtP->P_190) {
                B_5_0_0_idx_3 = _rtP->P_190;
            }
        }
        _rtB->B_27_205_0 = B_5_0_0_idx_3 + _rtDW->UnitDelay_DSTATE_b4;
          {
                 _rtB->B_27_206_0 = (_rtP->P_193[0])*_rtDW->DiscreteStateSpace_DSTATE_b[0] + (_rtP->P_193[1])*_rtDW->DiscreteStateSpace_DSTATE_b[1];
    _rtB->B_27_206_0 += _rtP->P_194*_rtB->B_27_205_0;
  }

   
  
        _rtB->B_27_208_0 = _rtP->P_196 * rtb_B_97_82_0;
        B_27_209_0 = ((B_27_79_0_idx_0 * rtb_B_27_80_0 + B_27_79_0_idx_1 * rtb_B_97_176_0) + B_27_79_0_idx_2 * rtb_B_97_396_0) * 0.66666666666666663;
        tmpForInput[0] = B_27_79_0_idx_0;
        tmpForInput[1] = B_27_79_0_idx_1;
        tmpForInput[2] = B_27_79_0_idx_2;
        tmpForInput[3] = rtb_B_27_80_0;
        tmpForInput[4] = rtb_B_27_81_0;
        tmpForInput[5] = rtb_B_97_176_0;
        tmpForInput[6] = rtb_B_97_226_0;
        tmpForInput[7] = rtb_B_97_396_0;
        tmpForInput[8] = rtb_B_97_177_0;
        rtb_B_97_177_0 = -0.0;
        for (i = 0; i < 9; i++) {
            rtb_B_97_177_0 += tmpForInput[i];
        }
        B_27_211_0 = _rtP->P_197 * rtb_B_97_177_0;
        tmpForInput[0] = B_27_79_0_idx_0;
        tmpForInput[1] = B_27_79_0_idx_1;
        tmpForInput[2] = B_27_79_0_idx_2;
        tmpForInput[3] = rtb_B_27_80_0;
        tmpForInput[4] = rtb_B_27_81_0;
        tmpForInput[5] = rtb_B_27_86_0;
        tmpForInput[6] = rtb_B_27_89_0;
        tmpForInput[7] = rtb_B_27_90_0;
        tmpForInput[8] = rtb_B_27_91_0;
        B_27_79_0_idx_0 = -0.0;
        for (i = 0; i < 9; i++) {
            B_27_79_0_idx_0 += tmpForInput[i];
        }
        B_27_213_0 = _rtP->P_198 * B_27_79_0_idx_0;
        tmpForInput[0] = B_27_95_0_idx_0;
        tmpForInput[1] = B_27_95_0_idx_1;
        tmpForInput[2] = B_27_95_0_idx_2;
        tmpForInput[3] = rtb_B_27_80_0;
        tmpForInput[4] = rtb_B_27_81_0;
        tmpForInput[5] = rtb_B_27_100_0;
        tmpForInput[6] = rtb_B_27_103_0;
        tmpForInput[7] = rtb_B_27_104_0;
        tmpForInput[8] = rtb_B_27_105_0;
        B_27_95_0_idx_0 = -0.0;
        for (i = 0; i < 9; i++) {
            B_27_95_0_idx_0 += tmpForInput[i];
        }
        B_27_215_0 = _rtP->P_199 * B_27_95_0_idx_0;
        rtb_B_97_82_0 = look1_pbinlxpw(muDoubleScalarRem(ssGetTaskTime(S, 2) + _rtP->P_200, _rtP->P_201) * _rtP->P_202, _rtP->P_204, _rtP->P_203, &_rtDW->m_bpIndex_b, 2U);
        B_27_224_0 = rtb_B_97_82_0 - _rtP->P_205;
    }
    i = ssIsSampleHit(S, 1, 0);
    if (i != 0) {
        _rtB->B_27_226_0 = _rtP->P_206;
    }
    rtb_B_97_82_0 = ssGetT(S) + _rtB->B_27_226_0;
    i = ssIsSampleHit(S, 1, 0);
    if (i != 0) {
        _rtB->B_27_228_0 = _rtP->P_207;
    }
    rtb_B_97_82_0 = look1_pbinlxpw(_rtP->P_208 * muDoubleScalarRem(rtb_B_97_82_0, _rtB->B_27_228_0), _rtP->P_210, _rtP->P_209, &_rtDW->m_bpIndex_j, 2U);
    i = ssIsSampleHit(S, 1, 0);
    if (i != 0) {
        _rtB->B_27_232_0 = _rtP->P_211;
    }
    B_27_233_0 = rtb_B_97_82_0 - _rtB->B_27_232_0;
    i = ssIsSampleHit(S, 4, 0);
    if (i != 0) {
        _rtB->B_27_234_0 = _rtDW->UnitDelay_DSTATE_f3;
    }
    _rtB->B_27_237_0 = _rtP->P_213 * _rtB->B_27_56_0;
    _rtB->B_27_238_0 = _rtP->P_214 * _rtB->B_27_50_0;
    _rtB->B_27_239_0 = _rtP->P_215 * _rtB->B_27_40_0;
    if (ssGetT(S) >= _rtP->P_217) {
        _rtB->B_27_241_0 = _rtB->B_27_47_0;
    } else {
        _rtB->B_27_241_0 = _rtB->B_27_240_0;
    }
    if (ssGetT(S) >= _rtP->P_219) {
        _rtB->B_27_243_0 = _rtB->B_27_47_0;
    } else {
        _rtB->B_27_243_0 = _rtB->B_27_242_0;
    }
    if (ssIsMajorTimeStep(S) != 0) {
        srUpdateBC(_rtDW->GridSupportingasCurrentSourceGridFeeding_SubsysRanBC);
    }
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    rtb_B_97_409_0 = _rtP->P_832 * _rtB->B_97_116_0[35];
    rtb_B_97_410_0 = _rtP->P_833 * _rtB->B_97_116_0[36];
    rtb_B_97_411_0 = _rtP->P_834 * _rtB->B_97_116_0[37];
    B_97_412_0_idx_0 = _rtP->P_835 * rtb_B_97_409_0;
    B_97_412_0_idx_1 = _rtP->P_835 * rtb_B_97_410_0;
    B_97_412_0_idx_2 = _rtP->P_835 * rtb_B_97_411_0;
}
rtPrevAction = _rtDW->If1_ActiveSubsystem;
if (ssIsMajorTimeStep(S) != 0) {
    rtAction = (int8_T)!(_rtB->B_97_405_0 > 0.5);
    _rtDW->If1_ActiveSubsystem = rtAction;
} else {
    rtAction = _rtDW->If1_ActiveSubsystem;
}
if (rtPrevAction != rtAction) {
    switch (rtPrevAction) {
      case 0:


      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_IfActionSubsystem2_Disable(S);      

        break;
      case 1:


      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_IfActionSubsystem2_Disable(S);      

        break;
    }
}
switch (rtAction) {
  case 0:
    if (rtAction != rtPrevAction) {


      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_IfActionSubsystem2_Enable(S);      

    }



      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_IfActionSubsystem2(S, _rtB->B_97_414_0, &B_97_418_0);
    if (ssIsMajorTimeStep(S) != 0) {
        srUpdateBC(_rtDW->IfActionSubsystem2.IfActionSubsystem2_SubsysRanBC);
    }
    break;
  case 1:
    if (rtAction != rtPrevAction) {


      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_IfActionSubsystem2_Enable(S);      

    }



      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_IfActionSubsystem2(S, _rtB->B_97_413_0, &B_97_418_0);
    if (ssIsMajorTimeStep(S) != 0) {
        srUpdateBC(_rtDW->IfActionSubsystem3.IfActionSubsystem2_SubsysRanBC);
    }
    break;
}
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    _rtB->B_97_419_0 = B_97_418_0;
    if (ssIsMajorTimeStep(S) != 0) {
        if (_rtB->B_97_419_0 > 0.0) {
            if (!_rtDW->Gridforming_MODE) {
                if (ssGetTaskTime(S, 1) != ssGetTStart(S)) {
                    ssSetBlockStateForSolverChangedAtMajorStep(S);
                }
    


      
      (void) memset(&(((XDis_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetContStateDisabled(S))->TransferFcn5_CSTATE), 0,
4*sizeof(boolean_T));

                _rtDW->Integ4_SYSTEM_ENABLE_bq = 1U;
                _rtDW->Gridforming_MODE = true;
            }
        } else {
            if (_rtDW->Gridforming_MODE) {
                ssSetBlockStateForSolverChangedAtMajorStep(S);
    


      
      (void) memset(&(((XDis_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetContStateDisabled(S))->TransferFcn5_CSTATE), 1,
4*sizeof(boolean_T));

                if (_rtDW->Subsystem1_MODE) {
                    _rtDW->Subsystem1_MODE = false;
                }
                if (_rtDW->Subsystempi2delay_MODE) {
                    _rtDW->Subsystempi2delay_MODE = false;
                }
                if (_rtDW->AutomaticGainControl_MODE_ab) {
                    _rtDW->AutomaticGainControl_MODE_ab = false;
                }
                _rtDW->Gridforming_MODE = false;
            }
        }
    }
}
if (_rtDW->Gridforming_MODE) {
    _rtB->B_36_1_0 = _rtB->B_36_0_0 - _rtB->B_97_363_0;
    i = ssIsSampleHit(S, 2, 0);
    if (i != 0) {
        B_5_0_0_idx_3 = _rtP->P_257 * _rtB->B_36_1_0 + _rtDW->DiscreteTimeIntegrator_DSTATE_b;
        if (B_5_0_0_idx_3 > _rtP->P_262) {
            B_5_0_0_idx_3 = _rtP->P_262;
        } else {
            if (B_5_0_0_idx_3 < _rtP->P_263) {
                B_5_0_0_idx_3 = _rtP->P_263;
            }
        }
        _rtB->B_36_10_0 = (B_5_0_0_idx_3 + _rtB->B_36_0_0) * ssGetTaskTime(S, 2) * _rtB->B_36_9_0;
    }
    _rtB->B_36_13_0 = 0.0;
    _rtB->B_36_13_0 += _rtP->P_268 * _rtX->TransferFcn5_CSTATE;
    _rtB->B_36_14_0 = _rtB->B_36_12_0 + _rtB->B_36_13_0;
    _rtB->B_36_15_0 = _rtB->B_36_10_0 + _rtB->B_36_14_0;
    _rtB->B_36_17_0 = muDoubleScalarMod(_rtB->B_36_15_0, _rtB->B_36_16_0);
    _rtB->B_36_19_0 = _rtB->B_36_18_0 - _rtB->B_97_362_0;
    i = ssIsSampleHit(S, 2, 0);
    if (i != 0) {
        B_5_0_0_idx_3 = _rtP->P_271 * _rtB->B_36_19_0 + _rtDW->DiscreteTimeIntegrator_DSTATE_m;
        if (B_5_0_0_idx_3 > _rtP->P_276) {
            B_5_0_0_idx_3 = _rtP->P_276;
        } else {
            if (B_5_0_0_idx_3 < _rtP->P_277) {
                B_5_0_0_idx_3 = _rtP->P_277;
            }
        }
        _rtB->B_36_26_0 = B_5_0_0_idx_3 + _rtB->B_36_25_0;
    }
    _rtB->B_36_27_0 = 0.0;
    _rtB->B_36_27_0 += _rtP->P_280 * _rtX->TransferFcn4_CSTATE_f;
    _rtB->B_36_28_0 = 0.0;
    _rtB->B_36_28_0 += _rtP->P_282 * _rtX->TransferFcn6_CSTATE;
    i = ssIsSampleHit(S, 1, 0);
    if (i != 0) {
        _rtB->B_36_32_0 = _rtB->B_36_31_0;
        if (ssIsMajorTimeStep(S) != 0) {
            if (_rtB->B_36_32_0 > 0) {
                if (!_rtDW->Subsystem1_MODE) {
                    if (ssGetTaskTime(S, 1) != ssGetTStart(S)) {
                        ssSetBlockStateForSolverChangedAtMajorStep(S);
                    }
  
                    _rtDW->Subsystem1_MODE = true;
                }
            } else {
                if (_rtDW->Subsystem1_MODE) {
                    ssSetBlockStateForSolverChangedAtMajorStep(S);
  
                    _rtDW->Subsystem1_MODE = false;
                }
            }
        }
    }
    if (_rtDW->Subsystem1_MODE) {
        _rtB->B_35_0_0 = _rtB->B_36_27_0 * muDoubleScalarCos(_rtB->B_36_17_0) - _rtB->B_36_28_0 * muDoubleScalarSin(_rtB->B_36_17_0);
        _rtB->B_35_1_0 = _rtB->B_36_27_0 * muDoubleScalarSin(_rtB->B_36_17_0) + _rtB->B_36_28_0 * muDoubleScalarCos(_rtB->B_36_17_0);
        if (ssIsMajorTimeStep(S) != 0) {
            srUpdateBC(_rtDW->Subsystem1_SubsysRanBC);
        }
    }
    i = ssIsSampleHit(S, 1, 0);
    if (i != 0) {
        _rtB->B_36_36_0 = _rtB->B_36_35_0;
        if (ssIsMajorTimeStep(S) != 0) {
            if (_rtB->B_36_36_0 > 0) {
                if (!_rtDW->Subsystempi2delay_MODE) {
                    if (ssGetTaskTime(S, 1) != ssGetTStart(S)) {
                        ssSetBlockStateForSolverChangedAtMajorStep(S);
                    }
  
                    _rtDW->Subsystempi2delay_MODE = true;
                }
            } else {
                if (_rtDW->Subsystempi2delay_MODE) {
                    ssSetBlockStateForSolverChangedAtMajorStep(S);
  
                    _rtDW->Subsystempi2delay_MODE = false;
                }
            }
        }
    }
    if (_rtDW->Subsystempi2delay_MODE) {
        _rtB->B_34_0_0 = _rtB->B_36_27_0 * muDoubleScalarSin(_rtB->B_36_17_0) + _rtB->B_36_28_0 * muDoubleScalarCos(_rtB->B_36_17_0);
        _rtB->B_34_1_0 = -_rtB->B_36_27_0 * muDoubleScalarCos(_rtB->B_36_17_0) + _rtB->B_36_28_0 * muDoubleScalarSin(_rtB->B_36_17_0);
        if (ssIsMajorTimeStep(S) != 0) {
            srUpdateBC(_rtDW->Subsystempi2delay_SubsysRanBC);
        }
    }
    if (_rtB->B_36_31_0 != 0) {
        _rtB->B_36_38_0[0] = _rtB->B_35_0_0;
        _rtB->B_36_38_0[1] = _rtB->B_35_1_0;
    } else {
        _rtB->B_36_38_0[0] = _rtB->B_34_0_0;
        _rtB->B_36_38_0[1] = _rtB->B_34_1_0;
    }
    _rtB->B_36_39_0 = 0.0;
    _rtB->B_36_39_0 += _rtP->P_287 * _rtX->TransferFcn7_CSTATE;
    _rtB->B_36_40_0[0] = _rtB->B_36_38_0[0];
    _rtB->B_36_40_0[1] = _rtB->B_36_38_0[1];
    _rtB->B_36_40_0[2] = _rtB->B_36_39_0;
    for (i = 0; i < 3; i++) {
        _rtB->B_36_41_0[i] = 0.0;
        _rtB->B_36_41_0[i] += _rtP->P_288[i] * _rtB->B_36_40_0[0];
        _rtB->B_36_41_0[i] += _rtP->P_288[i + 3] * _rtB->B_36_40_0[1];
        _rtB->B_36_41_0[i] += _rtP->P_288[i + 6] * _rtB->B_36_40_0[2];
    }
    _rtB->B_36_42_0[0] = _rtB->B_36_26_0 * _rtB->B_36_41_0[0];
    _rtB->B_36_42_0[1] = _rtB->B_36_26_0 * _rtB->B_36_41_0[1];
    _rtB->B_36_42_0[2] = _rtB->B_36_26_0 * _rtB->B_36_41_0[2];
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 36, 43, SS_CALL_MDL_OUTPUTS);
   
  
    rtb_B_27_80_0 = ssGetT(S);
    i = ssIsSampleHit(S, 2, 0);
    if (i != 0) {
        _rtB->B_36_46_0 = _rtP->P_289 * _rtB->B_36_19_0;
        _rtB->B_36_47_0 = _rtP->P_290 * _rtB->B_36_1_0;
    }
    B_36_48_0 = muDoubleScalarSin(_rtB->B_36_17_0);
    B_36_49_0 = muDoubleScalarCos(_rtB->B_36_17_0);
    i = ssIsSampleHit(S, 2, 0);
    if (i != 0) {
        _rtB->B_36_52_0 = muDoubleScalarMod(_rtDW->DiscreteTimeIntegrator_DSTATE_d, _rtB->B_36_51_0);
        rtb_B_27_81_0 = muDoubleScalarSin(_rtB->B_36_52_0);
        rtb_B_27_86_0 = muDoubleScalarCos(_rtB->B_36_52_0);
        rtb_B_27_89_0 = (0.0 - rtb_B_27_81_0 * _rtB->B_36_55_0) - rtb_B_27_86_0 * _rtB->B_36_57_0;
        rtb_B_27_90_0 = rtb_B_27_81_0 * _rtB->B_36_57_0 - rtb_B_27_86_0 * _rtB->B_36_55_0;
        rtb_B_27_91_0 = (0.0 - rtb_B_27_89_0) - rtb_B_27_81_0;
        rtb_B_27_100_0 = (0.0 - rtb_B_27_90_0) - rtb_B_27_86_0;
        _rtB->B_36_65_0 = ((B_97_412_0_idx_0 * rtb_B_27_86_0 + B_97_412_0_idx_1 * rtb_B_27_90_0) + B_97_412_0_idx_2 * rtb_B_27_100_0) * 0.66666666666666663;
        if (_rtDW->Integ4_SYSTEM_ENABLE_bq != 0) {
            _rtB->B_36_66_0 = _rtDW->Integ4_DSTATE_ba;
        } else {
            _rtB->B_36_66_0 = _rtP->P_296 * _rtB->B_36_65_0 + _rtDW->Integ4_DSTATE_ba;
        }
        rtb_B_27_103_0 = _rtDW->UnitDelay_DSTATE_lw;
        if (_rtDW->UnitDelay_DSTATE_lw > _rtP->P_299) {
            B_5_0_0_idx_3 = _rtP->P_299;
        } else if (_rtDW->UnitDelay_DSTATE_lw < _rtP->P_300) {
            B_5_0_0_idx_3 = _rtP->P_300;
        } else {
            B_5_0_0_idx_3 = _rtDW->UnitDelay_DSTATE_lw;
        }
        rtb_B_27_104_0 = 1.0 / B_5_0_0_idx_3 / 1.0e-5;
        rtb_B_27_105_0 = muDoubleScalarCeil(rtb_B_27_104_0);
        _rtB->B_36_71_0 = _rtP->P_301 * rtb_B_27_105_0;
                
          /* Level2 S-Function Block: '<S646>/B_36_3' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 36, 72, SS_CALL_MDL_OUTPUTS);

   
  
        rtb_B_97_82_0 = ssGetTaskTime(S, 2);
        if (rtb_B_97_82_0 < _rtP->P_303) {
            rtb_B_97_82_0 = _rtP->P_304;
        } else {
            rtb_B_97_82_0 = _rtP->P_305;
        }
        if (rtb_B_97_82_0 >= _rtP->P_307) {
            rtb_B_97_82_0 = rtb_B_27_104_0 - rtb_B_27_105_0;
            rtb_B_27_104_0 = ((_rtB->B_36_65_0 - _rtDW->UnitDelay_DSTATE_aq) * rtb_B_97_82_0 * _rtP->P_251 + _rtB->B_36_65_0) * (rtb_B_97_82_0 / rtb_B_27_104_0) + (_rtB->B_36_66_0 - _rtB->B_36_72_0) * rtb_B_27_103_0;
        } else {
            rtb_B_27_104_0 = _rtP->P_306;
        }
        _rtB->B_36_79_0 = _rtB->B_36_78_0;
        i = ssIsSampleHit(S, 2, 0);
        if ((i != 0) && (ssIsMajorTimeStep(S) != 0)) {
            if (_rtB->B_36_79_0 > 0.0) {
                if (!_rtDW->AutomaticGainControl_MODE_ab) {
                    if (ssGetTaskTime(S, 2) != ssGetTStart(S)) {
                        ssSetBlockStateForSolverChangedAtMajorStep(S);
                    }
  
                    _rtDW->Integ4_SYSTEM_ENABLE_lpf = 1U;
                    _rtDW->Integ4_SYSTEM_ENABLE_fw = 1U;
                    _rtDW->AutomaticGainControl_MODE_ab = true;
                }
            } else {
                if (_rtDW->AutomaticGainControl_MODE_ab) {
                    ssSetBlockStateForSolverChangedAtMajorStep(S);
  
                    _rtDW->AutomaticGainControl_MODE_ab = false;
                }
            }
        }
        if (_rtDW->AutomaticGainControl_MODE_ab) {
            i = ssIsSampleHit(S, 2, 0);
            if (i != 0) {
                rtb_B_27_105_0 = (0.0 - rtb_B_27_81_0 * _rtB->B_30_0_0) - rtb_B_27_86_0 * _rtB->B_30_2_0;
                rtb_B_97_177_0 = rtb_B_27_81_0 * _rtB->B_30_2_0 - rtb_B_27_86_0 * _rtB->B_30_0_0;
                B_27_79_0_idx_0 = (0.0 - rtb_B_27_105_0) - rtb_B_27_81_0;
                B_27_95_0_idx_0 = (0.0 - rtb_B_97_177_0) - rtb_B_27_86_0;
                _rtB->B_30_10_0 = ((B_97_412_0_idx_0 * rtb_B_27_81_0 + B_97_412_0_idx_1 * rtb_B_27_105_0) + B_97_412_0_idx_2 * B_27_79_0_idx_0) * 0.66666666666666663;
                if (_rtDW->Integ4_SYSTEM_ENABLE_lpf != 0) {
                    _rtB->B_30_11_0 = _rtDW->Integ4_DSTATE_kj;
                } else {
                    _rtB->B_30_11_0 = _rtP->P_225 * _rtB->B_30_10_0 + _rtDW->Integ4_DSTATE_kj;
                }
                if (rtb_B_27_103_0 > _rtP->P_227) {
                    rtb_B_97_82_0 = _rtP->P_227;
                } else if (rtb_B_27_103_0 < _rtP->P_228) {
                    rtb_B_97_82_0 = _rtP->P_228;
                } else {
                    rtb_B_97_82_0 = rtb_B_27_103_0;
                }
                B_5_0_0_idx_3 = 1.0 / rtb_B_97_82_0 / 1.0e-5;
                rtb_B_97_176_0 = muDoubleScalarCeil(B_5_0_0_idx_3);
                _rtB->B_30_15_0 = _rtP->P_229 * rtb_B_97_176_0;
                
          /* Level2 S-Function Block: '<S641>/B_30_0' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 30, 16, SS_CALL_MDL_OUTPUTS);

   
  
                rtb_B_97_82_0 = ssGetTaskTime(S, 2);
                if (rtb_B_97_82_0 < _rtP->P_231) {
                    rtb_B_97_82_0 = _rtP->P_232;
                } else {
                    rtb_B_97_82_0 = _rtP->P_233;
                }
                if (rtb_B_97_82_0 >= _rtP->P_235) {
                    rtb_B_97_82_0 = B_5_0_0_idx_3 - rtb_B_97_176_0;
                    B_5_0_0_idx_3 = ((_rtB->B_30_10_0 - _rtDW->UnitDelay_DSTATE_fa) * rtb_B_97_82_0 * _rtP->P_220 + _rtB->B_30_10_0) * (rtb_B_97_82_0 / B_5_0_0_idx_3) + (_rtB->B_30_11_0 - _rtB->B_30_16_0) * rtb_B_27_103_0;
                } else {
                    B_5_0_0_idx_3 = _rtP->P_234;
                }
                _rtB->B_30_22_0 = ((B_97_412_0_idx_0 * rtb_B_27_86_0 + B_97_412_0_idx_1 * rtb_B_97_177_0) + B_97_412_0_idx_2 * B_27_95_0_idx_0) * 0.66666666666666663;
                if (_rtDW->Integ4_SYSTEM_ENABLE_fw != 0) {
                    _rtB->B_30_23_0 = _rtDW->Integ4_DSTATE_d2;
                } else {
                    _rtB->B_30_23_0 = _rtP->P_236 * _rtB->B_30_22_0 + _rtDW->Integ4_DSTATE_d2;
                }
                if (rtb_B_27_103_0 > _rtP->P_238) {
                    rtb_B_97_82_0 = _rtP->P_238;
                } else if (rtb_B_27_103_0 < _rtP->P_239) {
                    rtb_B_97_82_0 = _rtP->P_239;
                } else {
                    rtb_B_97_82_0 = rtb_B_27_103_0;
                }
                rtb_B_97_176_0 = 1.0 / rtb_B_97_82_0 / 1.0e-5;
                rtb_B_97_226_0 = muDoubleScalarCeil(rtb_B_97_176_0);
                _rtB->B_30_27_0 = _rtP->P_240 * rtb_B_97_226_0;
                
          /* Level2 S-Function Block: '<S643>/B_30_2' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 30, 28, SS_CALL_MDL_OUTPUTS);

   
  
                rtb_B_97_82_0 = ssGetTaskTime(S, 2);
                if (rtb_B_97_82_0 < _rtP->P_242) {
                    rtb_B_97_82_0 = _rtP->P_243;
                } else {
                    rtb_B_97_82_0 = _rtP->P_244;
                }
                if (rtb_B_97_82_0 >= _rtP->P_246) {
                    rtb_B_97_82_0 = rtb_B_97_176_0 - rtb_B_97_226_0;
                    rtb_B_97_82_0 = ((_rtB->B_30_22_0 - _rtDW->UnitDelay_DSTATE_lu) * rtb_B_97_82_0 * _rtP->P_221 + _rtB->B_30_22_0) * (rtb_B_97_82_0 / rtb_B_97_176_0) + (_rtB->B_30_23_0 - _rtB->B_30_28_0) * rtb_B_27_103_0;
                } else {
                    rtb_B_97_82_0 = _rtP->P_245;
                }
                B_30_36_0 = _rtP->P_247 * muDoubleScalarAtan2(rtb_B_97_82_0, B_5_0_0_idx_3);
                tmpForInput[0] = B_97_412_0_idx_0;
                tmpForInput[1] = B_97_412_0_idx_1;
                tmpForInput[2] = B_97_412_0_idx_2;
                tmpForInput[3] = rtb_B_27_81_0;
                tmpForInput[4] = rtb_B_27_86_0;
                tmpForInput[5] = rtb_B_27_105_0;
                tmpForInput[6] = rtb_B_97_177_0;
                tmpForInput[7] = B_27_79_0_idx_0;
                tmpForInput[8] = B_27_95_0_idx_0;
                B_27_95_0_idx_0 = -0.0;
                for (i = 0; i < 9; i++) {
                    B_27_95_0_idx_0 += tmpForInput[i];
                }
                B_30_38_0 = _rtP->P_248 * B_27_95_0_idx_0;
                B_5_0_0_idx_3 = muDoubleScalarHypot(B_5_0_0_idx_3, rtb_B_97_82_0);
                if (B_5_0_0_idx_3 > _rtP->P_249) {
                    B_5_0_0_idx_3 = _rtP->P_249;
                } else {
                    if (B_5_0_0_idx_3 < _rtP->P_250) {
                        B_5_0_0_idx_3 = _rtP->P_250;
                    }
                }
                _rtB->B_30_40_0 = 1.0 / B_5_0_0_idx_3;
            }
            if (ssIsMajorTimeStep(S) != 0) {
                srUpdateBC(_rtDW->AutomaticGainControl_SubsysRanBC_p);
            }
        }
        rtb_B_97_82_0 = rtb_B_27_104_0 * _rtB->B_30_40_0;
        _rtB->B_36_86_0 = _rtP->P_314 * rtb_B_97_82_0 * _rtP->P_315;
        B_36_87_0 = _rtDW->UD_DSTATE_e;
        B_36_88_0 = _rtB->B_36_86_0 - B_36_87_0;
        B_5_0_0_idx_3 = (_rtP->P_309 * rtb_B_97_82_0 + _rtDW->DiscreteTimeIntegrator_DSTATE_ov) + B_36_88_0;
        if (B_5_0_0_idx_3 > _rtP->P_317) {
            _rtB->B_36_90_0 = _rtP->P_317;
        } else if (B_5_0_0_idx_3 < _rtP->P_318) {
            _rtB->B_36_90_0 = _rtP->P_318;
        } else {
            _rtB->B_36_90_0 = B_5_0_0_idx_3;
        }
        B_5_0_0_idx_3 = _rtP->P_319 * _rtB->B_36_90_0 - _rtDW->UnitDelay_DSTATE_n;
        if (B_5_0_0_idx_3 > _rtP->P_321) {
            B_5_0_0_idx_3 = _rtP->P_321;
        } else {
            if (B_5_0_0_idx_3 < _rtP->P_322) {
                B_5_0_0_idx_3 = _rtP->P_322;
            }
        }
        _rtB->B_36_95_0 = B_5_0_0_idx_3 + _rtDW->UnitDelay_DSTATE_n;
          {
                 _rtB->B_36_96_0 = (_rtP->P_325[0])*_rtDW->DiscreteStateSpace_DSTATE_pm[0] + (_rtP->P_325[1])*_rtDW->DiscreteStateSpace_DSTATE_pm[1];
    _rtB->B_36_96_0 += _rtP->P_326*_rtB->B_36_95_0;
  }

   
  
        _rtB->B_36_98_0 = _rtP->P_328 * rtb_B_97_82_0;
        B_36_99_0 = ((B_97_412_0_idx_0 * rtb_B_27_81_0 + B_97_412_0_idx_1 * rtb_B_27_89_0) + B_97_412_0_idx_2 * rtb_B_27_91_0) * 0.66666666666666663;
        tmpForInput[0] = B_97_412_0_idx_0;
        tmpForInput[1] = B_97_412_0_idx_1;
        tmpForInput[2] = B_97_412_0_idx_2;
        tmpForInput[3] = rtb_B_27_81_0;
        tmpForInput[4] = rtb_B_27_86_0;
        tmpForInput[5] = rtb_B_27_89_0;
        tmpForInput[6] = rtb_B_27_90_0;
        tmpForInput[7] = rtb_B_27_91_0;
        tmpForInput[8] = rtb_B_27_100_0;
        B_27_95_0_idx_1 = -0.0;
        for (i = 0; i < 9; i++) {
            B_27_95_0_idx_1 += tmpForInput[i];
        }
        B_36_101_0 = _rtP->P_329 * B_27_95_0_idx_1;
    }
    if (ssGetT(S) >= _rtP->P_330) {
        _rtB->B_36_102_0 = _rtB->B_36_13_0;
    } else {
        _rtB->B_36_102_0 = _rtB->B_36_52_0;
    }
    i = ssIsSampleHit(S, 2, 0);
    if (i != 0) {
        for (i = 0; i < 3; i++) {
            _rtB->B_36_104_0[i] = _rtP->P_332 * (_rtP->P_331[i + 6] * B_97_412_0_idx_2 + (_rtP->P_331[i + 3] * B_97_412_0_idx_1 + _rtP->P_331[i] * B_97_412_0_idx_0));
        }
        _rtB->B_36_108_0 = _rtB->B_36_107_0;
        if (_rtB->B_36_108_0 > 0) {
            _rtB->B_33_0_0 = _rtB->B_36_104_0[0] * muDoubleScalarCos(_rtB->B_36_52_0) + _rtB->B_36_104_0[1] * muDoubleScalarSin(_rtB->B_36_52_0);
            _rtB->B_33_1_0 = -_rtB->B_36_104_0[0] * muDoubleScalarSin(_rtB->B_36_52_0) + _rtB->B_36_104_0[1] * muDoubleScalarCos(_rtB->B_36_52_0);
            if (ssIsMajorTimeStep(S) != 0) {
                srUpdateBC(_rtDW->Subsystem1_SubsysRanBC_o);
            }
        }
        _rtB->B_36_112_0 = _rtB->B_36_111_0;
        if (_rtB->B_36_112_0 > 0) {
            _rtB->B_32_0_0 = _rtB->B_36_104_0[0] * muDoubleScalarSin(_rtB->B_36_52_0) - _rtB->B_36_104_0[1] * muDoubleScalarCos(_rtB->B_36_52_0);
            _rtB->B_32_1_0 = _rtB->B_36_104_0[0] * muDoubleScalarCos(_rtB->B_36_52_0) + _rtB->B_36_104_0[1] * muDoubleScalarSin(_rtB->B_36_52_0);
            if (ssIsMajorTimeStep(S) != 0) {
                srUpdateBC(_rtDW->Subsystempi2delay_SubsysRanBC_l);
            }
        }
        if (_rtB->B_36_107_0 != 0) {
            _rtB->B_36_114_0[0] = _rtB->B_33_0_0;
            _rtB->B_36_114_0[1] = _rtB->B_33_1_0;
        } else {
            _rtB->B_36_114_0[0] = _rtB->B_32_0_0;
            _rtB->B_36_114_0[1] = _rtB->B_32_1_0;
        }
    }
    if (rtb_B_27_80_0 >= _rtP->P_336) {
        _rtB->B_36_115_0[0] = _rtB->B_36_27_0;
        _rtB->B_36_115_0[1] = _rtB->B_36_28_0;
        _rtB->B_36_115_0[2] = _rtB->B_36_39_0;
    } else {
        _rtB->B_36_115_0[0] = _rtB->B_36_114_0[0];
        _rtB->B_36_115_0[1] = _rtB->B_36_114_0[1];
        _rtB->B_36_115_0[2] = _rtB->B_36_104_0[2];
    }
    if (ssIsMajorTimeStep(S) != 0) {
        srUpdateBC(_rtDW->Gridforming_SubsysRanBC);
    }
}
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 421, SS_CALL_MDL_OUTPUTS);
   
  
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 422, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 423, SS_CALL_MDL_OUTPUTS);
   
  
}
_rtB->B_97_425_0 = _rtP->P_839 * _rtB->B_97_361_0[0];
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    if (ssIsMajorTimeStep(S) != 0) {
        _rtDW->RelationalOperator_Mode[0] = (_rtB->B_97_361_0[1] > _rtB->B_97_425_0);
        _rtDW->RelationalOperator_Mode[1] = (_rtB->B_97_361_0[2] > _rtB->B_97_425_0);
    }
    _rtB->B_97_426_0[0] = _rtDW->RelationalOperator_Mode[0];
    _rtB->B_97_426_0[1] = _rtDW->RelationalOperator_Mode[1];
}
if (_rtB->B_97_426_0[0] >= _rtP->P_840) {
    B_97_427_0[0] = rtb_B_97_30_0[1];
} else {
    B_97_427_0[0] = _rtB->B_97_424_0;
}
if (_rtB->B_97_426_0[1] >= _rtP->P_840) {
    B_97_427_0[1] = rtb_B_97_30_0[2];
} else {
    B_97_427_0[1] = _rtB->B_97_424_0;
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_428_0 = _rtP->P_841 * _rtB->B_97_265_0[0];
}
rtb_B_97_82_0 = muDoubleScalarSin(_rtP->P_844 * ssGetTaskTime(S, 0) + _rtP->P_845) * _rtP->P_842 + _rtP->P_843;
_rtB->B_97_430_0 = rtb_B_97_82_0 * _rtB->B_97_428_0;
B_27_95_0_idx_1 = muDoubleScalarSin(_rtP->P_848 * ssGetTaskTime(S, 0) + _rtP->P_849) * _rtP->P_846 + _rtP->P_847;
_rtB->B_97_432_0 = _rtB->B_97_428_0 * B_27_95_0_idx_1;
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_433_0 = _rtP->P_850 * _rtB->B_97_265_0[1];
}
_rtB->B_97_434_0 = rtb_B_97_82_0 * _rtB->B_97_433_0;
_rtB->B_97_435_0 = _rtB->B_97_433_0 * B_27_95_0_idx_1;
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_436_0 = _rtP->P_851 * _rtB->B_97_265_0[2];
}
_rtB->B_97_437_0 = rtb_B_97_82_0 * _rtB->B_97_436_0;
_rtB->B_97_438_0 = _rtB->B_97_436_0 * B_27_95_0_idx_1;
B_97_439_0 = _rtP->P_852 * rtb_B_97_30_0[0];
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    B_97_291_1_idx_2 = _rtP->P_853 * _rtB->B_97_265_0[0];
    B_97_328_1_idx_1 = _rtP->P_853 * _rtB->B_97_265_0[1];
    B_97_328_1_idx_2 = _rtP->P_853 * _rtB->B_97_265_0[2];
    rtb_B_97_82_0 = muDoubleScalarMod(_rtDW->DiscreteTimeIntegrator_DSTATE_e, _rtB->B_97_442_0);
    B_27_95_0_idx_1 = muDoubleScalarSin(rtb_B_97_82_0);
    B_27_95_0_idx_0 = muDoubleScalarCos(rtb_B_97_82_0);
    rtb_B_27_80_0 = (0.0 - B_27_95_0_idx_1 * _rtB->B_97_446_0) - B_27_95_0_idx_0 * _rtB->B_97_448_0;
    rtb_B_27_81_0 = B_27_95_0_idx_1 * _rtB->B_97_448_0 - B_27_95_0_idx_0 * _rtB->B_97_446_0;
    rtb_B_27_86_0 = (0.0 - rtb_B_27_80_0) - B_27_95_0_idx_1;
    rtb_B_27_89_0 = (0.0 - rtb_B_27_81_0) - B_27_95_0_idx_0;
    _rtB->B_97_456_0 = ((B_97_291_1_idx_2 * B_27_95_0_idx_0 + B_97_328_1_idx_1 * rtb_B_27_81_0) + B_97_328_1_idx_2 * rtb_B_27_89_0) * 0.66666666666666663;
    if (_rtDW->Integ4_SYSTEM_ENABLE_e != 0) {
        _rtB->B_97_457_0 = _rtDW->Integ4_DSTATE_p;
    } else {
        _rtB->B_97_457_0 = _rtP->P_859 * _rtB->B_97_456_0 + _rtDW->Integ4_DSTATE_p;
    }
    _rtB->B_97_458_0 = _rtDW->UnitDelay_DSTATE_i;
    if (_rtB->B_97_458_0 > _rtP->P_862) {
        B_5_0_0_idx_3 = _rtP->P_862;
    } else if (_rtB->B_97_458_0 < _rtP->P_863) {
        B_5_0_0_idx_3 = _rtP->P_863;
    } else {
        B_5_0_0_idx_3 = _rtB->B_97_458_0;
    }
    rtb_B_27_90_0 = 1.0 / B_5_0_0_idx_3 / 1.0e-5;
    rtb_B_27_91_0 = muDoubleScalarCeil(rtb_B_27_90_0);
    _rtB->B_97_462_0 = _rtP->P_864 * rtb_B_27_91_0;
                
          /* Level2 S-Function Block: '<S686>/B_82_72' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 463, SS_CALL_MDL_OUTPUTS);

   
  
    rtb_B_97_82_0 = ssGetTaskTime(S, 2);
    if (rtb_B_97_82_0 < _rtP->P_866) {
        rtb_B_97_82_0 = _rtP->P_867;
    } else {
        rtb_B_97_82_0 = _rtP->P_868;
    }
    if (rtb_B_97_82_0 >= _rtP->P_870) {
        rtb_B_97_82_0 = rtb_B_27_90_0 - rtb_B_27_91_0;
        rtb_B_27_90_0 = ((_rtB->B_97_456_0 - _rtDW->UnitDelay_DSTATE_j) * rtb_B_97_82_0 * _rtP->P_368 + _rtB->B_97_456_0) * (rtb_B_97_82_0 / rtb_B_27_90_0) + (_rtB->B_97_457_0 - _rtB->B_97_463_0) * _rtB->B_97_458_0;
    } else {
        rtb_B_27_90_0 = _rtP->P_869;
    }
    _rtB->B_97_470_0 = _rtB->B_97_469_0;
    i = ssIsSampleHit(S, 2, 0);
    if ((i != 0) && (ssIsMajorTimeStep(S) != 0)) {
        if (_rtB->B_97_470_0 > 0.0) {
            if (!_rtDW->AutomaticGainControl_MODE_f) {
                if (ssGetTaskTime(S, 2) != ssGetTStart(S)) {
                    ssSetBlockStateForSolverChangedAtMajorStep(S);
                }
  
                _rtDW->Integ4_SYSTEM_ENABLE_ht = 1U;
                _rtDW->Integ4_SYSTEM_ENABLE_en = 1U;
                _rtDW->AutomaticGainControl_MODE_f = true;
            }
        } else {
            if (_rtDW->AutomaticGainControl_MODE_f) {
                ssSetBlockStateForSolverChangedAtMajorStep(S);
  
                _rtDW->AutomaticGainControl_MODE_f = false;
            }
        }
    }
    if (_rtDW->AutomaticGainControl_MODE_f) {
        i = ssIsSampleHit(S, 2, 0);
        if (i != 0) {
            rtb_B_27_91_0 = (0.0 - B_27_95_0_idx_1 * _rtB->B_44_0_0) - B_27_95_0_idx_0 * _rtB->B_44_2_0;
            rtb_B_27_100_0 = B_27_95_0_idx_1 * _rtB->B_44_2_0 - B_27_95_0_idx_0 * _rtB->B_44_0_0;
            rtb_B_27_103_0 = (0.0 - rtb_B_27_91_0) - B_27_95_0_idx_1;
            rtb_B_27_104_0 = (0.0 - rtb_B_27_100_0) - B_27_95_0_idx_0;
            _rtB->B_44_10_0 = ((B_97_291_1_idx_2 * B_27_95_0_idx_1 + B_97_328_1_idx_1 * rtb_B_27_91_0) + B_97_328_1_idx_2 * rtb_B_27_103_0) * 0.66666666666666663;
            if (_rtDW->Integ4_SYSTEM_ENABLE_ht != 0) {
                _rtB->B_44_11_0 = _rtDW->Integ4_DSTATE_ed;
            } else {
                _rtB->B_44_11_0 = _rtP->P_342 * _rtB->B_44_10_0 + _rtDW->Integ4_DSTATE_ed;
            }
            if (_rtB->B_97_458_0 > _rtP->P_344) {
                B_5_0_0_idx_3 = _rtP->P_344;
            } else if (_rtB->B_97_458_0 < _rtP->P_345) {
                B_5_0_0_idx_3 = _rtP->P_345;
            } else {
                B_5_0_0_idx_3 = _rtB->B_97_458_0;
            }
            rtb_B_27_105_0 = 1.0 / B_5_0_0_idx_3 / 1.0e-5;
            rtb_B_97_177_0 = muDoubleScalarCeil(rtb_B_27_105_0);
            _rtB->B_44_15_0 = _rtP->P_346 * rtb_B_97_177_0;
                
          /* Level2 S-Function Block: '<S681>/B_42_0' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 44, 16, SS_CALL_MDL_OUTPUTS);

   
  
            rtb_B_97_82_0 = ssGetTaskTime(S, 2);
            if (rtb_B_97_82_0 < _rtP->P_348) {
                rtb_B_97_82_0 = _rtP->P_349;
            } else {
                rtb_B_97_82_0 = _rtP->P_350;
            }
            if (rtb_B_97_82_0 >= _rtP->P_352) {
                rtb_B_97_82_0 = rtb_B_27_105_0 - rtb_B_97_177_0;
                rtb_B_27_105_0 = ((_rtB->B_44_10_0 - _rtDW->UnitDelay_DSTATE_g0) * rtb_B_97_82_0 * _rtP->P_337 + _rtB->B_44_10_0) * (rtb_B_97_82_0 / rtb_B_27_105_0) + (_rtB->B_44_11_0 - _rtB->B_44_16_0) * _rtB->B_97_458_0;
            } else {
                rtb_B_27_105_0 = _rtP->P_351;
            }
            _rtB->B_44_22_0 = ((B_97_291_1_idx_2 * B_27_95_0_idx_0 + B_97_328_1_idx_1 * rtb_B_27_100_0) + B_97_328_1_idx_2 * rtb_B_27_104_0) * 0.66666666666666663;
            if (_rtDW->Integ4_SYSTEM_ENABLE_en != 0) {
                _rtB->B_44_23_0 = _rtDW->Integ4_DSTATE_ok;
            } else {
                _rtB->B_44_23_0 = _rtP->P_353 * _rtB->B_44_22_0 + _rtDW->Integ4_DSTATE_ok;
            }
            if (_rtB->B_97_458_0 > _rtP->P_355) {
                B_5_0_0_idx_3 = _rtP->P_355;
            } else if (_rtB->B_97_458_0 < _rtP->P_356) {
                B_5_0_0_idx_3 = _rtP->P_356;
            } else {
                B_5_0_0_idx_3 = _rtB->B_97_458_0;
            }
            rtb_B_97_177_0 = 1.0 / B_5_0_0_idx_3 / 1.0e-5;
            B_27_79_0_idx_0 = muDoubleScalarCeil(rtb_B_97_177_0);
            _rtB->B_44_27_0 = _rtP->P_357 * B_27_79_0_idx_0;
                
          /* Level2 S-Function Block: '<S683>/B_42_2' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 44, 28, SS_CALL_MDL_OUTPUTS);

   
  
            rtb_B_97_82_0 = ssGetTaskTime(S, 2);
            if (rtb_B_97_82_0 < _rtP->P_359) {
                rtb_B_97_82_0 = _rtP->P_360;
            } else {
                rtb_B_97_82_0 = _rtP->P_361;
            }
            if (rtb_B_97_82_0 >= _rtP->P_363) {
                rtb_B_97_82_0 = rtb_B_97_177_0 - B_27_79_0_idx_0;
                rtb_B_97_82_0 = ((_rtB->B_44_22_0 - _rtDW->UnitDelay_DSTATE_il) * rtb_B_97_82_0 * _rtP->P_338 + _rtB->B_44_22_0) * (rtb_B_97_82_0 / rtb_B_97_177_0) + (_rtB->B_44_23_0 - _rtB->B_44_28_0) * _rtB->B_97_458_0;
            } else {
                rtb_B_97_82_0 = _rtP->P_362;
            }
            B_44_36_0 = _rtP->P_364 * muDoubleScalarAtan2(rtb_B_97_82_0, rtb_B_27_105_0);
            tmpForInput[0] = B_97_291_1_idx_2;
            tmpForInput[1] = B_97_328_1_idx_1;
            tmpForInput[2] = B_97_328_1_idx_2;
            tmpForInput[3] = B_27_95_0_idx_1;
            tmpForInput[4] = B_27_95_0_idx_0;
            tmpForInput[5] = rtb_B_27_91_0;
            tmpForInput[6] = rtb_B_27_100_0;
            tmpForInput[7] = rtb_B_27_103_0;
            tmpForInput[8] = rtb_B_27_104_0;
            rtb_B_27_91_0 = -0.0;
            for (i = 0; i < 9; i++) {
                rtb_B_27_91_0 += tmpForInput[i];
            }
            B_44_38_0 = _rtP->P_365 * rtb_B_27_91_0;
            B_5_0_0_idx_3 = muDoubleScalarHypot(rtb_B_27_105_0, rtb_B_97_82_0);
            if (B_5_0_0_idx_3 > _rtP->P_366) {
                B_5_0_0_idx_3 = _rtP->P_366;
            } else {
                if (B_5_0_0_idx_3 < _rtP->P_367) {
                    B_5_0_0_idx_3 = _rtP->P_367;
                }
            }
            _rtB->B_44_40_0 = 1.0 / B_5_0_0_idx_3;
        }
        if (ssIsMajorTimeStep(S) != 0) {
            srUpdateBC(_rtDW->AutomaticGainControl_SubsysRanBC_hn);
        }
    }
    rtb_B_97_82_0 = rtb_B_27_90_0 * _rtB->B_44_40_0;
    _rtB->B_97_477_0 = _rtP->P_877 * rtb_B_97_82_0 * _rtP->P_878;
    B_97_478_0 = _rtDW->UD_DSTATE;
    B_97_479_0 = _rtB->B_97_477_0 - B_97_478_0;
    B_5_0_0_idx_3 = (_rtP->P_872 * rtb_B_97_82_0 + _rtDW->DiscreteTimeIntegrator_DSTATE_j) + B_97_479_0;
    if (B_5_0_0_idx_3 > _rtP->P_880) {
        _rtB->B_97_481_0 = _rtP->P_880;
    } else if (B_5_0_0_idx_3 < _rtP->P_881) {
        _rtB->B_97_481_0 = _rtP->P_881;
    } else {
        _rtB->B_97_481_0 = B_5_0_0_idx_3;
    }
    B_5_0_0_idx_3 = _rtP->P_882 * _rtB->B_97_481_0 - _rtDW->UnitDelay_DSTATE_pc;
    if (B_5_0_0_idx_3 > _rtP->P_884) {
        B_5_0_0_idx_3 = _rtP->P_884;
    } else {
        if (B_5_0_0_idx_3 < _rtP->P_885) {
            B_5_0_0_idx_3 = _rtP->P_885;
        }
    }
    _rtB->B_97_486_0 = B_5_0_0_idx_3 + _rtDW->UnitDelay_DSTATE_pc;
          {
                 _rtB->B_97_487_0 = (_rtP->P_888[0])*_rtDW->DiscreteStateSpace_DSTATE[0] + (_rtP->P_888[1])*_rtDW->DiscreteStateSpace_DSTATE[1];
    _rtB->B_97_487_0 += _rtP->P_889*_rtB->B_97_486_0;
  }

   
  
    _rtB->B_97_489_0 = _rtP->P_891 * rtb_B_97_82_0;
    B_97_490_0 = ((B_97_291_1_idx_2 * B_27_95_0_idx_1 + B_97_328_1_idx_1 * rtb_B_27_80_0) + B_97_328_1_idx_2 * rtb_B_27_86_0) * 0.66666666666666663;
    tmpForInput[0] = B_97_291_1_idx_2;
    tmpForInput[1] = B_97_328_1_idx_1;
    tmpForInput[2] = B_97_328_1_idx_2;
    tmpForInput[3] = B_27_95_0_idx_1;
    tmpForInput[4] = B_27_95_0_idx_0;
    tmpForInput[5] = rtb_B_27_80_0;
    tmpForInput[6] = rtb_B_27_81_0;
    tmpForInput[7] = rtb_B_27_86_0;
    tmpForInput[8] = rtb_B_27_89_0;
    rtb_B_27_80_0 = -0.0;
    for (i = 0; i < 9; i++) {
        rtb_B_27_80_0 += tmpForInput[i];
    }
    B_97_492_0 = _rtP->P_892 * rtb_B_27_80_0;
    B_97_496_0[0] = !(_rtP->P_893 * _rtB->B_97_375_0[0] >= _rtB->B_97_382_0);
    B_97_496_0[1] = !(_rtP->P_893 * _rtB->B_97_375_0[1] >= _rtB->B_97_382_0);
    B_97_496_0[2] = !(_rtP->P_893 * _rtB->B_97_375_0[2] >= _rtB->B_97_382_0);
    rtb_B_97_82_0 = _rtP->P_896 * ssGetTaskTime(S, 2);
    rtb_B_27_80_0 = _rtP->P_898 * _rtP->P_897;
    B_97_291_1_idx_2 = muDoubleScalarSin((rtb_B_97_82_0 + rtb_B_27_80_0) + _rtP->P_894[0]) * _rtP->P_899;
    B_97_328_1_idx_1 = muDoubleScalarSin((rtb_B_97_82_0 + rtb_B_27_80_0) + _rtP->P_894[1]) * _rtP->P_899;
    B_97_328_1_idx_2 = muDoubleScalarSin((rtb_B_97_82_0 + rtb_B_27_80_0) + _rtP->P_894[2]) * _rtP->P_899;
    if (_rtP->P_895 >= _rtP->P_900) {
        B_97_507_0[0] = B_97_291_1_idx_2;
        B_97_507_0[1] = B_97_328_1_idx_1;
        B_97_507_0[2] = B_97_328_1_idx_2;
    } else {
        B_97_507_0[0] = 0.0;
        B_97_507_0[1] = 0.0;
        B_97_507_0[2] = 0.0;
    }
}
if (_rtB->B_97_406_0 >= _rtP->P_901) {
    _rtB->B_97_508_0[0] = _rtB->B_27_151_0[0];
    _rtB->B_97_508_0[1] = _rtB->B_27_151_0[1];
    _rtB->B_97_508_0[2] = _rtB->B_27_151_0[2];
} else {
    _rtB->B_97_508_0[0] = _rtB->B_36_42_0[0];
    _rtB->B_97_508_0[1] = _rtB->B_36_42_0[1];
    _rtB->B_97_508_0[2] = _rtB->B_36_42_0[2];
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    B_97_509_0[0] = _rtP->P_902 * rtb_B_97_409_0;
    B_97_509_0[1] = _rtP->P_902 * rtb_B_97_410_0;
    B_97_509_0[2] = _rtP->P_902 * rtb_B_97_411_0;
    B_97_291_1_idx_2 = _rtP->P_903 * _rtB->B_97_116_0[29] * _rtP->P_906;
    B_97_328_1_idx_1 = _rtP->P_904 * _rtB->B_97_116_0[30] * _rtP->P_906;
    B_97_328_1_idx_2 = _rtP->P_905 * _rtB->B_97_116_0[31] * _rtP->P_906;
    rtb_B_97_82_0 = _rtP->P_907 * ssGetTaskTime(S, 2);
    rtb_B_97_574_0 = muDoubleScalarSin(rtb_B_97_82_0);
    rtb_B_97_575_0 = muDoubleScalarCos(rtb_B_97_82_0);
    rtb_B_97_580_0 = (0.0 - rtb_B_97_574_0 * _rtB->B_97_576_0) - rtb_B_97_575_0 * _rtB->B_97_578_0;
    rtb_B_97_583_0 = rtb_B_97_574_0 * _rtB->B_97_578_0 - rtb_B_97_575_0 * _rtB->B_97_576_0;
    rtb_B_97_584_0 = (0.0 - rtb_B_97_580_0) - rtb_B_97_574_0;
    rtb_B_97_585_0 = (0.0 - rtb_B_97_583_0) - rtb_B_97_575_0;
    _rtB->B_97_586_0 = ((B_97_291_1_idx_2 * rtb_B_97_574_0 + B_97_328_1_idx_1 * rtb_B_97_580_0) + B_97_328_1_idx_2 * rtb_B_97_584_0) * 0.66666666666666663;
    if (_rtDW->Integ4_SYSTEM_ENABLE_j != 0) {
        _rtB->B_97_587_0 = _rtDW->Integ4_DSTATE_lr;
    } else {
        _rtB->B_97_587_0 = _rtP->P_910 * _rtB->B_97_586_0 + _rtDW->Integ4_DSTATE_lr;
    }
    _rtB->B_97_588_0 = _rtP->P_912;
                
          /* Level2 S-Function Block: '<S182>/B_82_134' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 589, SS_CALL_MDL_OUTPUTS);

   
  
    rtb_B_97_82_0 = ssGetTaskTime(S, 2);
    if (rtb_B_97_82_0 < _rtP->P_915) {
        rtb_B_97_82_0 = _rtP->P_916;
    } else {
        rtb_B_97_82_0 = _rtP->P_917;
    }
    if (rtb_B_97_82_0 >= _rtP->P_919) {
        rtb_B_97_409_0 = (_rtB->B_97_587_0 - _rtB->B_97_589_0) * _rtP->P_913 + (_rtP->P_19 * _rtB->B_97_586_0 - _rtP->P_18 * _rtDW->UnitDelay_DSTATE_dl);
    } else {
        rtb_B_97_409_0 = _rtP->P_918;
    }
    _rtB->B_97_596_0 = ((B_97_291_1_idx_2 * rtb_B_97_575_0 + B_97_328_1_idx_1 * rtb_B_97_583_0) + B_97_328_1_idx_2 * rtb_B_97_585_0) * 0.66666666666666663;
    if (_rtDW->Integ4_SYSTEM_ENABLE_fv != 0) {
        _rtB->B_97_597_0 = _rtDW->Integ4_DSTATE_g;
    } else {
        _rtB->B_97_597_0 = _rtP->P_920 * _rtB->B_97_596_0 + _rtDW->Integ4_DSTATE_g;
    }
    _rtB->B_97_598_0 = _rtP->P_922;
                
          /* Level2 S-Function Block: '<S183>/B_82_136' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 599, SS_CALL_MDL_OUTPUTS);

   
  
    rtb_B_97_82_0 = ssGetTaskTime(S, 2);
    if (rtb_B_97_82_0 < _rtP->P_925) {
        rtb_B_97_82_0 = _rtP->P_926;
    } else {
        rtb_B_97_82_0 = _rtP->P_927;
    }
    if (rtb_B_97_82_0 >= _rtP->P_929) {
        rtb_B_97_410_0 = (_rtB->B_97_597_0 - _rtB->B_97_599_0) * _rtP->P_923 + (_rtP->P_21 * _rtB->B_97_596_0 - _rtP->P_20 * _rtDW->UnitDelay_DSTATE_e0);
    } else {
        rtb_B_97_410_0 = _rtP->P_928;
    }
    B_97_321_0_idx_0 = _rtP->P_930 * _rtB->B_97_116_0[51] * _rtP->P_933;
    B_97_321_0_idx_1 = _rtP->P_931 * _rtB->B_97_116_0[52] * _rtP->P_933;
    B_97_321_0_idx_2 = _rtP->P_932 * _rtB->B_97_116_0[53] * _rtP->P_933;
    rtb_B_97_82_0 = _rtP->P_934 * ssGetTaskTime(S, 2);
    rtb_B_97_613_0 = muDoubleScalarSin(rtb_B_97_82_0);
    rtb_B_97_614_0 = muDoubleScalarCos(rtb_B_97_82_0);
    rtb_B_97_619_0 = (0.0 - rtb_B_97_613_0 * _rtB->B_97_615_0) - rtb_B_97_614_0 * _rtB->B_97_617_0;
    rtb_B_97_622_0 = rtb_B_97_613_0 * _rtB->B_97_617_0 - rtb_B_97_614_0 * _rtB->B_97_615_0;
    rtb_B_97_623_0 = (0.0 - rtb_B_97_619_0) - rtb_B_97_613_0;
    rtb_B_97_624_0 = (0.0 - rtb_B_97_622_0) - rtb_B_97_614_0;
    _rtB->B_97_625_0 = ((B_97_321_0_idx_0 * rtb_B_97_613_0 + B_97_321_0_idx_1 * rtb_B_97_619_0) + B_97_321_0_idx_2 * rtb_B_97_623_0) * 0.66666666666666663;
    if (_rtDW->Integ4_SYSTEM_ENABLE_jr != 0) {
        _rtB->B_97_626_0 = _rtDW->Integ4_DSTATE_j;
    } else {
        _rtB->B_97_626_0 = _rtP->P_937 * _rtB->B_97_625_0 + _rtDW->Integ4_DSTATE_j;
    }
    _rtB->B_97_627_0 = _rtP->P_939;
                
          /* Level2 S-Function Block: '<S176>/B_82_138' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 628, SS_CALL_MDL_OUTPUTS);

   
  
    rtb_B_97_82_0 = ssGetTaskTime(S, 2);
    if (rtb_B_97_82_0 < _rtP->P_942) {
        rtb_B_97_82_0 = _rtP->P_943;
    } else {
        rtb_B_97_82_0 = _rtP->P_944;
    }
    if (rtb_B_97_82_0 >= _rtP->P_946) {
        rtb_B_97_411_0 = (_rtB->B_97_626_0 - _rtB->B_97_628_0) * _rtP->P_940 + (_rtP->P_15 * _rtB->B_97_625_0 - _rtP->P_14 * _rtDW->UnitDelay_DSTATE_k);
    } else {
        rtb_B_97_411_0 = _rtP->P_945;
    }
    _rtB->B_97_635_0 = ((B_97_321_0_idx_0 * rtb_B_97_614_0 + B_97_321_0_idx_1 * rtb_B_97_622_0) + B_97_321_0_idx_2 * rtb_B_97_624_0) * 0.66666666666666663;
    if (_rtDW->Integ4_SYSTEM_ENABLE_g != 0) {
        _rtB->B_97_636_0 = _rtDW->Integ4_DSTATE_b;
    } else {
        _rtB->B_97_636_0 = _rtP->P_947 * _rtB->B_97_635_0 + _rtDW->Integ4_DSTATE_b;
    }
    _rtB->B_97_637_0 = _rtP->P_949;
                
          /* Level2 S-Function Block: '<S177>/B_82_140' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 638, SS_CALL_MDL_OUTPUTS);

   
  
    rtb_B_97_82_0 = ssGetTaskTime(S, 2);
    if (rtb_B_97_82_0 < _rtP->P_952) {
        rtb_B_97_82_0 = _rtP->P_953;
    } else {
        rtb_B_97_82_0 = _rtP->P_954;
    }
    if (rtb_B_97_82_0 >= _rtP->P_956) {
        rtb_B_27_80_0 = (_rtB->B_97_636_0 - _rtB->B_97_638_0) * _rtP->P_950 + (_rtP->P_17 * _rtB->B_97_635_0 - _rtP->P_16 * _rtDW->UnitDelay_DSTATE_g);
    } else {
        rtb_B_27_80_0 = _rtP->P_955;
    }
    rtb_B_97_82_0 = muDoubleScalarHypot(rtb_B_97_409_0, rtb_B_97_410_0) * muDoubleScalarHypot(rtb_B_97_411_0, rtb_B_27_80_0);
    rtb_B_97_409_0 = (_rtP->P_957 * muDoubleScalarAtan2(rtb_B_97_410_0, rtb_B_97_409_0) - _rtP->P_958 * muDoubleScalarAtan2(rtb_B_27_80_0, rtb_B_97_411_0)) * _rtP->P_959;
    _rtB->B_97_657_0[0] = rtb_B_97_82_0 * muDoubleScalarCos(rtb_B_97_409_0) * _rtP->P_960;
    _rtB->B_97_657_0[1] = muDoubleScalarSin(rtb_B_97_409_0) * rtb_B_97_82_0 * _rtP->P_960;
    _rtB->B_97_658_0[0] = _rtP->P_961 * _rtB->B_97_657_0[0];
    _rtB->B_97_658_0[1] = _rtP->P_961 * _rtB->B_97_657_0[1];
    ssCallAccelRunBlock(S, 97, 659, SS_CALL_MDL_OUTPUTS);
    _rtB->B_97_295_0[0] = _rtP->P_962 * _rtB->B_97_116_0[45] * _rtP->P_965;
    _rtB->B_97_295_0[1] = _rtP->P_963 * _rtB->B_97_116_0[46] * _rtP->P_965;
    _rtB->B_97_295_0[2] = _rtP->P_964 * _rtB->B_97_116_0[47] * _rtP->P_965;
    if (_rtDW->systemEnable_k != 0) {
        _rtDW->lastSin_l = muDoubleScalarSin(_rtP->P_968 * ssGetTaskTime(S, 2));
        _rtDW->lastCos_n = muDoubleScalarCos(_rtP->P_968 * ssGetTaskTime(S, 2));
        _rtDW->systemEnable_k = 0;
    }
    rtb_B_97_82_0 = ((_rtDW->lastSin_l * _rtP->P_972 + _rtDW->lastCos_n * _rtP->P_971) * _rtP->P_970 + (_rtDW->lastCos_n * _rtP->P_972 - _rtDW->lastSin_l * _rtP->P_971) * _rtP->P_969) * _rtP->P_966 + _rtP->P_967;
    _rtB->B_97_665_0[0] = _rtB->B_97_295_0[0] * rtb_B_97_82_0;
    _rtB->B_97_665_0[1] = _rtB->B_97_295_0[1] * rtb_B_97_82_0;
    _rtB->B_97_665_0[2] = _rtB->B_97_295_0[2] * rtb_B_97_82_0;
    if (_rtDW->Integ4_SYSTEM_ENABLE_fl != 0) {
        _rtB->B_97_666_0[0] = _rtDW->Integ4_DSTATE_a[0];
        _rtB->B_97_666_0[1] = _rtDW->Integ4_DSTATE_a[1];
        _rtB->B_97_666_0[2] = _rtDW->Integ4_DSTATE_a[2];
    } else {
        _rtB->B_97_666_0[0] = _rtP->P_973 * _rtB->B_97_665_0[0] + _rtDW->Integ4_DSTATE_a[0];
        _rtB->B_97_666_0[1] = _rtP->P_973 * _rtB->B_97_665_0[1] + _rtDW->Integ4_DSTATE_a[1];
        _rtB->B_97_666_0[2] = _rtP->P_973 * _rtB->B_97_665_0[2] + _rtDW->Integ4_DSTATE_a[2];
    }
    _rtB->B_97_667_0 = _rtP->P_975;
                
          /* Level2 S-Function Block: '<S778>/B_82_142' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 668, SS_CALL_MDL_OUTPUTS);

   
  
    if (ssGetTaskTime(S, 2) >= _rtB->B_97_667_0) {
        _rtB->B_97_675_0[0] = (_rtB->B_97_666_0[0] - _rtB->B_97_668_0[0]) * _rtP->P_984 + (_rtP->P_452 * _rtB->B_97_665_0[0] - _rtP->P_451 * _rtDW->UnitDelay_DSTATE_g5[0]);
        _rtB->B_97_675_0[1] = (_rtB->B_97_666_0[1] - _rtB->B_97_668_0[1]) * _rtP->P_984 + (_rtP->P_452 * _rtB->B_97_665_0[1] - _rtP->P_451 * _rtDW->UnitDelay_DSTATE_g5[1]);
        _rtB->B_97_675_0[2] = (_rtB->B_97_666_0[2] - _rtB->B_97_668_0[2]) * _rtP->P_984 + (_rtP->P_452 * _rtB->B_97_665_0[2] - _rtP->P_451 * _rtDW->UnitDelay_DSTATE_g5[2]);
    } else {
        _rtB->B_97_675_0[0] = _rtDW->UnitDelay1_DSTATE_g[0];
        _rtB->B_97_675_0[1] = _rtDW->UnitDelay1_DSTATE_g[1];
        _rtB->B_97_675_0[2] = _rtDW->UnitDelay1_DSTATE_g[2];
    }
    if (_rtDW->systemEnable_o != 0) {
        _rtDW->lastSin_j = muDoubleScalarSin(_rtP->P_989 * ssGetTaskTime(S, 2));
        _rtDW->lastCos_ne = muDoubleScalarCos(_rtP->P_989 * ssGetTaskTime(S, 2));
        _rtDW->systemEnable_o = 0;
    }
    rtb_B_97_82_0 = ((_rtDW->lastSin_j * _rtP->P_993 + _rtDW->lastCos_ne * _rtP->P_992) * _rtP->P_991 + (_rtDW->lastCos_ne * _rtP->P_993 - _rtDW->lastSin_j * _rtP->P_992) * _rtP->P_990) * _rtP->P_987 + _rtP->P_988;
    _rtB->B_97_677_0[0] = _rtB->B_97_295_0[0] * rtb_B_97_82_0;
    _rtB->B_97_677_0[1] = _rtB->B_97_295_0[1] * rtb_B_97_82_0;
    _rtB->B_97_677_0[2] = _rtB->B_97_295_0[2] * rtb_B_97_82_0;
    if (_rtDW->Integ4_SYSTEM_ENABLE_ei != 0) {
        _rtB->B_97_678_0[0] = _rtDW->Integ4_DSTATE_e[0];
        _rtB->B_97_678_0[1] = _rtDW->Integ4_DSTATE_e[1];
        _rtB->B_97_678_0[2] = _rtDW->Integ4_DSTATE_e[2];
    } else {
        _rtB->B_97_678_0[0] = _rtP->P_994 * _rtB->B_97_677_0[0] + _rtDW->Integ4_DSTATE_e[0];
        _rtB->B_97_678_0[1] = _rtP->P_994 * _rtB->B_97_677_0[1] + _rtDW->Integ4_DSTATE_e[1];
        _rtB->B_97_678_0[2] = _rtP->P_994 * _rtB->B_97_677_0[2] + _rtDW->Integ4_DSTATE_e[2];
    }
    _rtB->B_97_679_0 = _rtP->P_996;
                
          /* Level2 S-Function Block: '<S776>/B_82_144' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 680, SS_CALL_MDL_OUTPUTS);

   
  
    if (ssGetTaskTime(S, 2) >= _rtB->B_97_679_0) {
        _rtB->B_97_687_0[0] = (_rtB->B_97_678_0[0] - _rtB->B_97_680_0[0]) * _rtP->P_1005 + (_rtP->P_450 * _rtB->B_97_677_0[0] - _rtP->P_449 * _rtDW->UnitDelay_DSTATE_b[0]);
        _rtB->B_97_687_0[1] = (_rtB->B_97_678_0[1] - _rtB->B_97_680_0[1]) * _rtP->P_1005 + (_rtP->P_450 * _rtB->B_97_677_0[1] - _rtP->P_449 * _rtDW->UnitDelay_DSTATE_b[1]);
        _rtB->B_97_687_0[2] = (_rtB->B_97_678_0[2] - _rtB->B_97_680_0[2]) * _rtP->P_1005 + (_rtP->P_450 * _rtB->B_97_677_0[2] - _rtP->P_449 * _rtDW->UnitDelay_DSTATE_b[2]);
    } else {
        _rtB->B_97_687_0[0] = _rtDW->UnitDelay1_DSTATE_j2[0];
        _rtB->B_97_687_0[1] = _rtDW->UnitDelay1_DSTATE_j2[1];
        _rtB->B_97_687_0[2] = _rtDW->UnitDelay1_DSTATE_j2[2];
    }
    rtb_B_97_66_0 = muDoubleScalarHypot(_rtB->B_97_675_0[0], _rtB->B_97_687_0[0]);
    rtb_B_97_30_0[0] = muDoubleScalarAtan2(_rtB->B_97_687_0[0], _rtB->B_97_675_0[0]);
    B_5_0_0_idx_3 = muDoubleScalarHypot(_rtB->B_97_675_0[1], _rtB->B_97_687_0[1]);
    rtb_B_97_30_0[1] = muDoubleScalarAtan2(_rtB->B_97_687_0[1], _rtB->B_97_675_0[1]);
    B_97_726_1_idx_0 = _rtB->B_97_675_0[2];
    B_97_719_0_idx_2 = _rtB->B_97_687_0[2];
    rtb_B_97_30_0[2] = muDoubleScalarAtan2(_rtB->B_97_687_0[2], _rtB->B_97_675_0[2]);
    _rtB->B_97_265_0[0] = _rtP->P_1008 * _rtB->B_97_116_0[65] * _rtP->P_1011;
    _rtB->B_97_265_0[1] = _rtP->P_1009 * _rtB->B_97_116_0[66] * _rtP->P_1011;
    _rtB->B_97_265_0[2] = _rtP->P_1010 * _rtB->B_97_116_0[67] * _rtP->P_1011;
    if (_rtDW->systemEnable_i != 0) {
        _rtDW->lastSin_h = muDoubleScalarSin(_rtP->P_1014 * ssGetTaskTime(S, 2));
        _rtDW->lastCos_n3 = muDoubleScalarCos(_rtP->P_1014 * ssGetTaskTime(S, 2));
        _rtDW->systemEnable_i = 0;
    }
    rtb_B_97_82_0 = ((_rtDW->lastSin_h * _rtP->P_1018 + _rtDW->lastCos_n3 * _rtP->P_1017) * _rtP->P_1016 + (_rtDW->lastCos_n3 * _rtP->P_1018 - _rtDW->lastSin_h * _rtP->P_1017) * _rtP->P_1015) * _rtP->P_1012 + _rtP->P_1013;
    _rtB->B_97_695_0[0] = _rtB->B_97_265_0[0] * rtb_B_97_82_0;
    _rtB->B_97_695_0[1] = _rtB->B_97_265_0[1] * rtb_B_97_82_0;
    _rtB->B_97_695_0[2] = _rtB->B_97_265_0[2] * rtb_B_97_82_0;
    if (_rtDW->Integ4_SYSTEM_ENABLE_f0 != 0) {
        _rtB->B_97_696_0[0] = _rtDW->Integ4_DSTATE_i[0];
        _rtB->B_97_696_0[1] = _rtDW->Integ4_DSTATE_i[1];
        _rtB->B_97_696_0[2] = _rtDW->Integ4_DSTATE_i[2];
    } else {
        _rtB->B_97_696_0[0] = _rtP->P_1019 * _rtB->B_97_695_0[0] + _rtDW->Integ4_DSTATE_i[0];
        _rtB->B_97_696_0[1] = _rtP->P_1019 * _rtB->B_97_695_0[1] + _rtDW->Integ4_DSTATE_i[1];
        _rtB->B_97_696_0[2] = _rtP->P_1019 * _rtB->B_97_695_0[2] + _rtDW->Integ4_DSTATE_i[2];
    }
    _rtB->B_97_697_0 = _rtP->P_1021;
                
          /* Level2 S-Function Block: '<S784>/B_82_146' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 698, SS_CALL_MDL_OUTPUTS);

   
  
    if (ssGetTaskTime(S, 2) >= _rtB->B_97_697_0) {
        _rtB->B_97_705_0[0] = (_rtB->B_97_696_0[0] - _rtB->B_97_698_0[0]) * _rtP->P_1030 + (_rtP->P_456 * _rtB->B_97_695_0[0] - _rtP->P_455 * _rtDW->UnitDelay_DSTATE_dlx[0]);
        _rtB->B_97_705_0[1] = (_rtB->B_97_696_0[1] - _rtB->B_97_698_0[1]) * _rtP->P_1030 + (_rtP->P_456 * _rtB->B_97_695_0[1] - _rtP->P_455 * _rtDW->UnitDelay_DSTATE_dlx[1]);
        _rtB->B_97_705_0[2] = (_rtB->B_97_696_0[2] - _rtB->B_97_698_0[2]) * _rtP->P_1030 + (_rtP->P_456 * _rtB->B_97_695_0[2] - _rtP->P_455 * _rtDW->UnitDelay_DSTATE_dlx[2]);
    } else {
        _rtB->B_97_705_0[0] = _rtDW->UnitDelay1_DSTATE_j5[0];
        _rtB->B_97_705_0[1] = _rtDW->UnitDelay1_DSTATE_j5[1];
        _rtB->B_97_705_0[2] = _rtDW->UnitDelay1_DSTATE_j5[2];
    }
    if (_rtDW->systemEnable_or != 0) {
        _rtDW->lastSin_n = muDoubleScalarSin(_rtP->P_1035 * ssGetTaskTime(S, 2));
        _rtDW->lastCos_d = muDoubleScalarCos(_rtP->P_1035 * ssGetTaskTime(S, 2));
        _rtDW->systemEnable_or = 0;
    }
    rtb_B_97_82_0 = ((_rtDW->lastSin_n * _rtP->P_1039 + _rtDW->lastCos_d * _rtP->P_1038) * _rtP->P_1037 + (_rtDW->lastCos_d * _rtP->P_1039 - _rtDW->lastSin_n * _rtP->P_1038) * _rtP->P_1036) * _rtP->P_1033 + _rtP->P_1034;
    _rtB->B_97_707_0[0] = _rtB->B_97_265_0[0] * rtb_B_97_82_0;
    _rtB->B_97_707_0[1] = _rtB->B_97_265_0[1] * rtb_B_97_82_0;
    _rtB->B_97_707_0[2] = _rtB->B_97_265_0[2] * rtb_B_97_82_0;
    if (_rtDW->Integ4_SYSTEM_ENABLE_m != 0) {
        _rtB->B_97_708_0[0] = _rtDW->Integ4_DSTATE_kd[0];
        _rtB->B_97_708_0[1] = _rtDW->Integ4_DSTATE_kd[1];
        _rtB->B_97_708_0[2] = _rtDW->Integ4_DSTATE_kd[2];
    } else {
        _rtB->B_97_708_0[0] = _rtP->P_1040 * _rtB->B_97_707_0[0] + _rtDW->Integ4_DSTATE_kd[0];
        _rtB->B_97_708_0[1] = _rtP->P_1040 * _rtB->B_97_707_0[1] + _rtDW->Integ4_DSTATE_kd[1];
        _rtB->B_97_708_0[2] = _rtP->P_1040 * _rtB->B_97_707_0[2] + _rtDW->Integ4_DSTATE_kd[2];
    }
    _rtB->B_97_709_0 = _rtP->P_1042;
                
          /* Level2 S-Function Block: '<S782>/B_82_148' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 710, SS_CALL_MDL_OUTPUTS);

   
  
    if (ssGetTaskTime(S, 2) >= _rtB->B_97_709_0) {
        _rtB->B_97_717_0[0] = (_rtB->B_97_708_0[0] - _rtB->B_97_710_0[0]) * _rtP->P_1051 + (_rtP->P_454 * _rtB->B_97_707_0[0] - _rtP->P_453 * _rtDW->UnitDelay_DSTATE_m[0]);
        _rtB->B_97_717_0[1] = (_rtB->B_97_708_0[1] - _rtB->B_97_710_0[1]) * _rtP->P_1051 + (_rtP->P_454 * _rtB->B_97_707_0[1] - _rtP->P_453 * _rtDW->UnitDelay_DSTATE_m[1]);
        _rtB->B_97_717_0[2] = (_rtB->B_97_708_0[2] - _rtB->B_97_710_0[2]) * _rtP->P_1051 + (_rtP->P_454 * _rtB->B_97_707_0[2] - _rtP->P_453 * _rtDW->UnitDelay_DSTATE_m[2]);
    } else {
        _rtB->B_97_717_0[0] = _rtDW->UnitDelay1_DSTATE_a[0];
        _rtB->B_97_717_0[1] = _rtDW->UnitDelay1_DSTATE_a[1];
        _rtB->B_97_717_0[2] = _rtDW->UnitDelay1_DSTATE_a[2];
    }
    B_97_719_0_idx_0 = rtb_B_97_66_0 * muDoubleScalarHypot(_rtB->B_97_705_0[0], _rtB->B_97_717_0[0]) * _rtP->P_1054;
    rtb_B_97_30_0[0] = (_rtP->P_1055 * rtb_B_97_30_0[0] - _rtP->P_1056 * muDoubleScalarAtan2(_rtB->B_97_717_0[0], _rtB->B_97_705_0[0])) * _rtP->P_1057;
    B_97_719_0_idx_1 = B_5_0_0_idx_3 * muDoubleScalarHypot(_rtB->B_97_705_0[1], _rtB->B_97_717_0[1]) * _rtP->P_1054;
    rtb_B_97_30_0[1] = (_rtP->P_1055 * rtb_B_97_30_0[1] - _rtP->P_1056 * muDoubleScalarAtan2(_rtB->B_97_717_0[1], _rtB->B_97_705_0[1])) * _rtP->P_1057;
    B_97_719_0_idx_2 = muDoubleScalarHypot(B_97_726_1_idx_0, B_97_719_0_idx_2) * muDoubleScalarHypot(_rtB->B_97_705_0[2], _rtB->B_97_717_0[2]) * _rtP->P_1054;
    rtb_B_97_30_0[2] = (_rtP->P_1055 * rtb_B_97_30_0[2] - _rtP->P_1056 * muDoubleScalarAtan2(_rtB->B_97_717_0[2], _rtB->B_97_705_0[2])) * _rtP->P_1057;
    muDoubleScalarSinCos(rtb_B_97_30_0[0], &B_97_726_0_idx_0, &B_97_726_1_idx_0);
    muDoubleScalarSinCos(rtb_B_97_30_0[1], &B_97_726_0_idx_1, &B_97_726_1_idx_1);
    muDoubleScalarSinCos(rtb_B_97_30_0[2], &B_97_726_0_idx_2, &B_97_726_1_idx_2);
    B_97_726_1_idx_0 *= B_97_719_0_idx_0;
    B_97_726_1_idx_1 *= B_97_719_0_idx_1;
    B_97_726_1_idx_2 *= B_97_719_0_idx_2;
    _rtB->B_97_728_0 = (B_97_726_1_idx_0 + B_97_726_1_idx_1) + B_97_726_1_idx_2;
}
_rtB->B_97_729_0 = _rtX->integ1_CSTATE_nd;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T_PWORK_bx.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T_PWORK_bx.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_1059;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        _rtB->B_97_730_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T_IWORK_a.CircularBufSize,
        &_rtDW->T_IWORK_a.Last,
        _rtDW->T_IWORK_a.Tail,
        _rtDW->T_IWORK_a.Head,
        _rtP->P_1060,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
_rtB->B_97_731_0 = _rtB->B_97_729_0 - _rtB->B_97_730_0;
_rtB->B_97_732_0 = _rtX->Integ2_CSTATE_i;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK_o.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK_o.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_1062;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        _rtB->B_97_733_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T1_IWORK_n.CircularBufSize,
        &_rtDW->T1_IWORK_n.Last,
        _rtDW->T1_IWORK_n.Tail,
        _rtDW->T1_IWORK_n.Head,
        _rtP->P_1063,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
_rtB->B_97_734_0 = _rtB->B_97_732_0 - _rtB->B_97_733_0;
_rtB->B_97_735_0.re = _rtB->B_97_731_0;
_rtB->B_97_735_0.im = _rtB->B_97_734_0;
_rtB->B_97_736_0 = _rtX->integ1_CSTATE_d;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T_PWORK_le.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T_PWORK_le.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_1065;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        _rtB->B_97_737_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T_IWORK_c.CircularBufSize,
        &_rtDW->T_IWORK_c.Last,
        _rtDW->T_IWORK_c.Tail,
        _rtDW->T_IWORK_c.Head,
        _rtP->P_1066,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
_rtB->B_97_738_0 = _rtB->B_97_736_0 - _rtB->B_97_737_0;
_rtB->B_97_739_0 = _rtX->Integ2_CSTATE_pa;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK_c.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK_c.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_1068;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        _rtB->B_97_740_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T1_IWORK_c.CircularBufSize,
        &_rtDW->T1_IWORK_c.Last,
        _rtDW->T1_IWORK_c.Tail,
        _rtDW->T1_IWORK_c.Head,
        _rtP->P_1069,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
_rtB->B_97_741_0 = _rtB->B_97_739_0 - _rtB->B_97_740_0;
_rtB->B_97_742_0.re = _rtB->B_97_738_0;
_rtB->B_97_742_0.im = _rtB->B_97_741_0;
_rtB->B_97_743_0 = _rtX->integ1_CSTATE_e;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T_PWORK_i.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T_PWORK_i.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_1071;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        _rtB->B_97_744_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T_IWORK_o.CircularBufSize,
        &_rtDW->T_IWORK_o.Last,
        _rtDW->T_IWORK_o.Tail,
        _rtDW->T_IWORK_o.Head,
        _rtP->P_1072,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
_rtB->B_97_745_0 = _rtB->B_97_743_0 - _rtB->B_97_744_0;
_rtB->B_97_746_0 = _rtX->Integ2_CSTATE_e;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK_h.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK_h.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_1074;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        _rtB->B_97_747_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T1_IWORK_j.CircularBufSize,
        &_rtDW->T1_IWORK_j.Last,
        _rtDW->T1_IWORK_j.Tail,
        _rtDW->T1_IWORK_j.Head,
        _rtP->P_1075,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
_rtB->B_97_748_0 = _rtB->B_97_746_0 - _rtB->B_97_747_0;
_rtB->B_97_749_0.re = _rtB->B_97_745_0;
_rtB->B_97_749_0.im = _rtB->B_97_748_0;
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    _rtB->B_97_751_0 = _rtB->B_97_750_0;
}


      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_NegSeqComputation(S, _rtB->B_97_751_0, _rtB->B_97_735_0, _rtB->B_97_742_0, _rtB->B_97_749_0, &_rtB->PosSeqComputation_f, &_rtDW->PosSeqComputation_f, &_rtP->PosSeqComputation_f);
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    _rtB->B_97_754_0 = _rtB->B_97_753_0;
}


      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_NegSeqComputation(S, _rtB->B_97_754_0, _rtB->B_97_735_0, _rtB->B_97_742_0, _rtB->B_97_749_0, &_rtB->NegSeqComputation_b, &_rtDW->NegSeqComputation_b, &_rtP->NegSeqComputation_b);
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    _rtB->B_97_757_0 = _rtB->B_97_756_0;
}


      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_ZeroSeqComputation(S, _rtB->B_97_757_0, _rtB->B_97_735_0, _rtB->B_97_742_0, _rtB->B_97_749_0, &_rtB->ZeroSeqComputation_i, &_rtDW->ZeroSeqComputation_i, &_rtP->ZeroSeqComputation_i);
_rtB->B_97_759_0[0] = muDoubleScalarHypot(_rtB->PosSeqComputation_f.B_39_2_0.re, _rtB->PosSeqComputation_f.B_39_2_0.im);
_rtB->B_97_759_0[1] = muDoubleScalarHypot(_rtB->NegSeqComputation_b.B_39_2_0.re, _rtB->NegSeqComputation_b.B_39_2_0.im);
_rtB->B_97_759_0[2] = muDoubleScalarHypot(_rtB->ZeroSeqComputation_i.B_41_1_0.re, _rtB->ZeroSeqComputation_i.B_41_1_0.im);
rtb_B_97_30_0[0] = muDoubleScalarAtan2(_rtB->PosSeqComputation_f.B_39_2_0.im, _rtB->PosSeqComputation_f.B_39_2_0.re);
rtb_B_97_30_0[1] = muDoubleScalarAtan2(_rtB->NegSeqComputation_b.B_39_2_0.im, _rtB->NegSeqComputation_b.B_39_2_0.re);
rtb_B_97_30_0[2] = muDoubleScalarAtan2(_rtB->ZeroSeqComputation_i.B_41_1_0.im, _rtB->ZeroSeqComputation_i.B_41_1_0.re);
_rtB->B_97_760_0 = _rtP->P_1079 * _rtB->B_97_759_0[0];
_rtB->B_97_761_0 = 0.0;
_rtB->B_97_761_0 += _rtP->P_1081 * _rtX->TransferFcn1_CSTATE_k;
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    B_97_719_0_idx_0 *= B_97_726_0_idx_0;
    B_97_719_0_idx_1 *= B_97_726_0_idx_1;
    _rtB->B_97_763_0 = (B_97_719_0_idx_0 + B_97_719_0_idx_1) + B_97_719_0_idx_2 * B_97_726_0_idx_2;
}
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 764, SS_CALL_MDL_OUTPUTS);
   
  
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    for (i = 0; i < 6; i++) {
        if (_rtB->B_97_116_1[i + 19] >= _rtP->P_1084) {
            B_5_0_0_idx_3 = _rtB->B_97_116_0[i + 19] * _rtP->P_1083;
        } else {
            B_5_0_0_idx_3 = _rtP->P_1082;
        }
        if (B_5_0_0_idx_3 > _rtP->P_1085) {
            _rtB->B_97_768_0[i] = _rtP->P_1085;
        } else if (B_5_0_0_idx_3 < _rtP->P_1086) {
            _rtB->B_97_768_0[i] = _rtP->P_1086;
        } else {
            _rtB->B_97_768_0[i] = B_5_0_0_idx_3;
        }
    }
    for (i = 0; i < 6; i++) {
        if (_rtDW->UnitDelay_DSTATE_c[i] != 0.0) {
            B_97_772_0[i] = _rtB->B_97_770_0[i];
        } else {
            B_97_772_0[i] = _rtB->B_97_771_0[i];
        }
    }
    rtb_B_97_82_0 = _rtP->P_1091 * ssGetTaskTime(S, 2);
    rtb_B_97_777_0 = _rtP->P_1093 * _rtP->P_1092;
    B_97_726_1_idx_0 = _rtP->P_1094[0];
    B_97_781_0_idx_0 = muDoubleScalarSin((rtb_B_97_82_0 + rtb_B_97_777_0) + _rtP->P_1094[0]) * _rtP->P_1090;
    B_97_726_1_idx_1 = _rtP->P_1094[1];
    B_97_781_0_idx_1 = muDoubleScalarSin((rtb_B_97_82_0 + rtb_B_97_777_0) + _rtP->P_1094[1]) * _rtP->P_1090;
    B_97_726_1_idx_2 = _rtP->P_1094[2];
    B_97_781_0_idx_2 = muDoubleScalarSin((rtb_B_97_82_0 + rtb_B_97_777_0) + _rtP->P_1094[2]) * _rtP->P_1090;
    rtb_B_97_777_0 = _rtP->P_1095;
}
i = ssIsSampleHit(S, 3, 0);
if (i != 0) {
    _rtB->B_97_783_0[0] = _rtDW->UnitDelay4_DSTATE[0];
    _rtB->B_97_783_0[1] = _rtDW->UnitDelay4_DSTATE[1];
    _rtB->B_97_783_0[2] = _rtDW->UnitDelay4_DSTATE[2];
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    if (rtb_B_97_777_0 >= _rtP->P_1097) {
        B_97_726_1_idx_0 = B_97_781_0_idx_0;
        B_97_726_1_idx_1 = B_97_781_0_idx_1;
        B_97_726_1_idx_2 = B_97_781_0_idx_2;
    } else {
        B_97_726_1_idx_0 = _rtB->B_97_783_0[0];
        B_97_726_1_idx_1 = _rtB->B_97_783_0[1];
        B_97_726_1_idx_2 = _rtB->B_97_783_0[2];
    }
}
_rtB->B_97_791_0 = look1_binlxpw(muDoubleScalarRem(ssGetT(S) - _rtB->B_97_786_0, _rtB->B_97_788_0), _rtP->P_1100, _rtP->P_1099, 3U);
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    rtb_B_97_384_0 = (B_97_726_1_idx_0 >= _rtB->B_97_791_0);
    rtb_B_97_76_0 = !rtb_B_97_384_0;
    rtb_B_97_111_0 = rtb_B_97_384_0;
    rtb_B_97_384_0 = (B_97_726_1_idx_1 >= _rtB->B_97_791_0);
    B_97_385_0_idx_1 = !rtb_B_97_384_0;
    B_97_384_0_idx_1 = rtb_B_97_384_0;
    rtb_B_97_384_0 = (B_97_726_1_idx_2 >= _rtB->B_97_791_0);
    _rtB->B_97_796_0[0] = rtb_B_97_111_0;
    _rtB->B_97_796_0[1] = rtb_B_97_76_0;
    _rtB->B_97_796_0[2] = B_97_384_0_idx_1;
    _rtB->B_97_796_0[3] = B_97_385_0_idx_1;
    _rtB->B_97_796_0[4] = rtb_B_97_384_0;
    _rtB->B_97_796_0[5] = !rtb_B_97_384_0;
    tmpForInput[0] = B_97_321_0_idx_0;
    tmpForInput[1] = B_97_321_0_idx_1;
    tmpForInput[2] = B_97_321_0_idx_2;
    tmpForInput[3] = rtb_B_97_613_0;
    tmpForInput[4] = rtb_B_97_614_0;
    tmpForInput[5] = rtb_B_97_619_0;
    tmpForInput[6] = rtb_B_97_622_0;
    tmpForInput[7] = rtb_B_97_623_0;
    tmpForInput[8] = rtb_B_97_624_0;
    rtb_B_97_803_0 = -0.0;
    for (i = 0; i < 9; i++) {
        rtb_B_97_803_0 += tmpForInput[i];
    }
    B_97_804_0 = _rtP->P_1101 * rtb_B_97_803_0;
    tmpForInput[0] = B_97_291_1_idx_2;
    tmpForInput[1] = B_97_328_1_idx_1;
    tmpForInput[2] = B_97_328_1_idx_2;
    tmpForInput[3] = rtb_B_97_574_0;
    tmpForInput[4] = rtb_B_97_575_0;
    tmpForInput[5] = rtb_B_97_580_0;
    tmpForInput[6] = rtb_B_97_583_0;
    tmpForInput[7] = rtb_B_97_584_0;
    tmpForInput[8] = rtb_B_97_585_0;
    rtb_B_97_805_0 = -0.0;
    for (i = 0; i < 9; i++) {
        rtb_B_97_805_0 += tmpForInput[i];
    }
    B_97_806_0 = _rtP->P_1102 * rtb_B_97_805_0;
}
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 807, SS_CALL_MDL_OUTPUTS);
   
  
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 808, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 809, SS_CALL_MDL_OUTPUTS);
   
  
}
_rtB->B_97_811_0 = _rtP->P_1104 * _rtB->B_97_759_0[0];
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    if (ssIsMajorTimeStep(S) != 0) {
        _rtDW->RelationalOperator_Mode_e[0] = (_rtB->B_97_759_0[1] > _rtB->B_97_811_0);
        _rtDW->RelationalOperator_Mode_e[1] = (_rtB->B_97_759_0[2] > _rtB->B_97_811_0);
    }
    _rtB->B_97_812_0[0] = _rtDW->RelationalOperator_Mode_e[0];
    _rtB->B_97_812_0[1] = _rtDW->RelationalOperator_Mode_e[1];
}
if (_rtB->B_97_812_0[0] >= _rtP->P_1105) {
    B_97_813_0[0] = rtb_B_97_30_0[1];
} else {
    B_97_813_0[0] = _rtB->B_97_810_0;
}
if (_rtB->B_97_812_0[1] >= _rtP->P_1105) {
    B_97_813_0[1] = rtb_B_97_30_0[2];
} else {
    B_97_813_0[1] = _rtB->B_97_810_0;
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_814_0 = _rtP->P_1106 * _rtB->B_97_295_0[0];
}
rtb_B_97_82_0 = muDoubleScalarSin(_rtP->P_1109 * ssGetTaskTime(S, 0) + _rtP->P_1110) * _rtP->P_1107 + _rtP->P_1108;
_rtB->B_97_816_0 = rtb_B_97_82_0 * _rtB->B_97_814_0;
rtb_B_97_574_0 = muDoubleScalarSin(_rtP->P_1113 * ssGetTaskTime(S, 0) + _rtP->P_1114) * _rtP->P_1111 + _rtP->P_1112;
_rtB->B_97_818_0 = _rtB->B_97_814_0 * rtb_B_97_574_0;
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_819_0 = _rtP->P_1115 * _rtB->B_97_295_0[1];
}
_rtB->B_97_820_0 = rtb_B_97_82_0 * _rtB->B_97_819_0;
_rtB->B_97_821_0 = _rtB->B_97_819_0 * rtb_B_97_574_0;
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_822_0 = _rtP->P_1116 * _rtB->B_97_295_0[2];
}
_rtB->B_97_823_0 = rtb_B_97_82_0 * _rtB->B_97_822_0;
_rtB->B_97_824_0 = _rtB->B_97_822_0 * rtb_B_97_574_0;
B_97_825_0 = _rtP->P_1117 * rtb_B_97_30_0[0];
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    B_97_719_0_idx_0 = _rtP->P_1118 * _rtB->B_97_295_0[0];
    B_97_719_0_idx_1 = _rtP->P_1118 * _rtB->B_97_295_0[1];
    B_97_719_0_idx_2 = _rtP->P_1118 * _rtB->B_97_295_0[2];
    rtb_B_97_82_0 = muDoubleScalarMod(_rtDW->DiscreteTimeIntegrator_DSTATE_p, _rtB->B_97_828_0);
    rtb_B_97_574_0 = muDoubleScalarSin(rtb_B_97_82_0);
    rtb_B_97_803_0 = muDoubleScalarCos(rtb_B_97_82_0);
    rtb_B_97_805_0 = (0.0 - rtb_B_97_574_0 * _rtB->B_97_832_0) - rtb_B_97_803_0 * _rtB->B_97_834_0;
    rtb_B_97_575_0 = rtb_B_97_574_0 * _rtB->B_97_834_0 - rtb_B_97_803_0 * _rtB->B_97_832_0;
    rtb_B_97_580_0 = (0.0 - rtb_B_97_805_0) - rtb_B_97_574_0;
    rtb_B_97_583_0 = (0.0 - rtb_B_97_575_0) - rtb_B_97_803_0;
    _rtB->B_97_842_0 = ((B_97_719_0_idx_0 * rtb_B_97_803_0 + B_97_719_0_idx_1 * rtb_B_97_575_0) + B_97_719_0_idx_2 * rtb_B_97_583_0) * 0.66666666666666663;
    if (_rtDW->Integ4_SYSTEM_ENABLE_i != 0) {
        _rtB->B_97_843_0 = _rtDW->Integ4_DSTATE_ay;
    } else {
        _rtB->B_97_843_0 = _rtP->P_1124 * _rtB->B_97_842_0 + _rtDW->Integ4_DSTATE_ay;
    }
    _rtB->B_97_844_0 = _rtDW->UnitDelay_DSTATE_l;
    if (_rtB->B_97_844_0 > _rtP->P_1127) {
        B_5_0_0_idx_3 = _rtP->P_1127;
    } else if (_rtB->B_97_844_0 < _rtP->P_1128) {
        B_5_0_0_idx_3 = _rtP->P_1128;
    } else {
        B_5_0_0_idx_3 = _rtB->B_97_844_0;
    }
    rtb_B_97_584_0 = 1.0 / B_5_0_0_idx_3 / 1.0e-5;
    rtb_B_97_585_0 = muDoubleScalarCeil(rtb_B_97_584_0);
    _rtB->B_97_848_0 = _rtP->P_1129 * rtb_B_97_585_0;
                
          /* Level2 S-Function Block: '<S770>/B_82_169' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 849, SS_CALL_MDL_OUTPUTS);

   
  
    rtb_B_97_82_0 = ssGetTaskTime(S, 2);
    if (rtb_B_97_82_0 < _rtP->P_1131) {
        rtb_B_97_82_0 = _rtP->P_1132;
    } else {
        rtb_B_97_82_0 = _rtP->P_1133;
    }
    if (rtb_B_97_82_0 >= _rtP->P_1135) {
        rtb_B_97_82_0 = rtb_B_97_584_0 - rtb_B_97_585_0;
        rtb_B_97_584_0 = ((_rtB->B_97_842_0 - _rtDW->UnitDelay_DSTATE_dw) * rtb_B_97_82_0 * _rtP->P_448 + _rtB->B_97_842_0) * (rtb_B_97_82_0 / rtb_B_97_584_0) + (_rtB->B_97_843_0 - _rtB->B_97_849_0) * _rtB->B_97_844_0;
    } else {
        rtb_B_97_584_0 = _rtP->P_1134;
    }
    _rtB->B_97_856_0 = _rtB->B_97_855_0;
    i = ssIsSampleHit(S, 2, 0);
    if ((i != 0) && (ssIsMajorTimeStep(S) != 0)) {
        if (_rtB->B_97_856_0 > 0.0) {
            if (!_rtDW->AutomaticGainControl_MODE_p) {
                if (ssGetTaskTime(S, 2) != ssGetTStart(S)) {
                    ssSetBlockStateForSolverChangedAtMajorStep(S);
                }
  
                _rtDW->Integ4_SYSTEM_ENABLE_k5 = 1U;
                _rtDW->Integ4_SYSTEM_ENABLE_ks = 1U;
                _rtDW->AutomaticGainControl_MODE_p = true;
            }
        } else {
            if (_rtDW->AutomaticGainControl_MODE_p) {
                ssSetBlockStateForSolverChangedAtMajorStep(S);
  
                _rtDW->AutomaticGainControl_MODE_p = false;
            }
        }
    }
    if (_rtDW->AutomaticGainControl_MODE_p) {
        i = ssIsSampleHit(S, 2, 0);
        if (i != 0) {
            rtb_B_97_585_0 = (0.0 - rtb_B_97_574_0 * _rtB->B_66_0_0) - rtb_B_97_803_0 * _rtB->B_66_2_0;
            rtb_B_97_613_0 = rtb_B_97_574_0 * _rtB->B_66_2_0 - rtb_B_97_803_0 * _rtB->B_66_0_0;
            rtb_B_97_614_0 = (0.0 - rtb_B_97_585_0) - rtb_B_97_574_0;
            rtb_B_97_619_0 = (0.0 - rtb_B_97_613_0) - rtb_B_97_803_0;
            _rtB->B_66_10_0 = ((B_97_719_0_idx_0 * rtb_B_97_574_0 + B_97_719_0_idx_1 * rtb_B_97_585_0) + B_97_719_0_idx_2 * rtb_B_97_614_0) * 0.66666666666666663;
            if (_rtDW->Integ4_SYSTEM_ENABLE_k5 != 0) {
                _rtB->B_66_11_0 = _rtDW->Integ4_DSTATE_o;
            } else {
                _rtB->B_66_11_0 = _rtP->P_422 * _rtB->B_66_10_0 + _rtDW->Integ4_DSTATE_o;
            }
            if (_rtB->B_97_844_0 > _rtP->P_424) {
                B_5_0_0_idx_3 = _rtP->P_424;
            } else if (_rtB->B_97_844_0 < _rtP->P_425) {
                B_5_0_0_idx_3 = _rtP->P_425;
            } else {
                B_5_0_0_idx_3 = _rtB->B_97_844_0;
            }
            rtb_B_97_622_0 = 1.0 / B_5_0_0_idx_3 / 1.0e-5;
            rtb_B_97_623_0 = muDoubleScalarCeil(rtb_B_97_622_0);
            _rtB->B_66_15_0 = _rtP->P_426 * rtb_B_97_623_0;
                
          /* Level2 S-Function Block: '<S765>/B_58_0' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 66, 16, SS_CALL_MDL_OUTPUTS);

   
  
            rtb_B_97_82_0 = ssGetTaskTime(S, 2);
            if (rtb_B_97_82_0 < _rtP->P_428) {
                rtb_B_97_82_0 = _rtP->P_429;
            } else {
                rtb_B_97_82_0 = _rtP->P_430;
            }
            if (rtb_B_97_82_0 >= _rtP->P_432) {
                rtb_B_97_82_0 = rtb_B_97_622_0 - rtb_B_97_623_0;
                rtb_B_97_622_0 = ((_rtB->B_66_10_0 - _rtDW->UnitDelay_DSTATE_a) * rtb_B_97_82_0 * _rtP->P_417 + _rtB->B_66_10_0) * (rtb_B_97_82_0 / rtb_B_97_622_0) + (_rtB->B_66_11_0 - _rtB->B_66_16_0) * _rtB->B_97_844_0;
            } else {
                rtb_B_97_622_0 = _rtP->P_431;
            }
            _rtB->B_66_22_0 = ((B_97_719_0_idx_0 * rtb_B_97_803_0 + B_97_719_0_idx_1 * rtb_B_97_613_0) + B_97_719_0_idx_2 * rtb_B_97_619_0) * 0.66666666666666663;
            if (_rtDW->Integ4_SYSTEM_ENABLE_ks != 0) {
                _rtB->B_66_23_0 = _rtDW->Integ4_DSTATE_g1;
            } else {
                _rtB->B_66_23_0 = _rtP->P_433 * _rtB->B_66_22_0 + _rtDW->Integ4_DSTATE_g1;
            }
            if (_rtB->B_97_844_0 > _rtP->P_435) {
                B_5_0_0_idx_3 = _rtP->P_435;
            } else if (_rtB->B_97_844_0 < _rtP->P_436) {
                B_5_0_0_idx_3 = _rtP->P_436;
            } else {
                B_5_0_0_idx_3 = _rtB->B_97_844_0;
            }
            rtb_B_97_623_0 = 1.0 / B_5_0_0_idx_3 / 1.0e-5;
            rtb_B_97_624_0 = muDoubleScalarCeil(rtb_B_97_623_0);
            _rtB->B_66_27_0 = _rtP->P_437 * rtb_B_97_624_0;
                
          /* Level2 S-Function Block: '<S767>/B_58_2' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 66, 28, SS_CALL_MDL_OUTPUTS);

   
  
            rtb_B_97_82_0 = ssGetTaskTime(S, 2);
            if (rtb_B_97_82_0 < _rtP->P_439) {
                rtb_B_97_82_0 = _rtP->P_440;
            } else {
                rtb_B_97_82_0 = _rtP->P_441;
            }
            if (rtb_B_97_82_0 >= _rtP->P_443) {
                rtb_B_97_82_0 = rtb_B_97_623_0 - rtb_B_97_624_0;
                rtb_B_97_82_0 = ((_rtB->B_66_22_0 - _rtDW->UnitDelay_DSTATE_b0) * rtb_B_97_82_0 * _rtP->P_418 + _rtB->B_66_22_0) * (rtb_B_97_82_0 / rtb_B_97_623_0) + (_rtB->B_66_23_0 - _rtB->B_66_28_0) * _rtB->B_97_844_0;
            } else {
                rtb_B_97_82_0 = _rtP->P_442;
            }
            B_66_36_0 = _rtP->P_444 * muDoubleScalarAtan2(rtb_B_97_82_0, rtb_B_97_622_0);
            tmpForInput[0] = B_97_719_0_idx_0;
            tmpForInput[1] = B_97_719_0_idx_1;
            tmpForInput[2] = B_97_719_0_idx_2;
            tmpForInput[3] = rtb_B_97_574_0;
            tmpForInput[4] = rtb_B_97_803_0;
            tmpForInput[5] = rtb_B_97_585_0;
            tmpForInput[6] = rtb_B_97_613_0;
            tmpForInput[7] = rtb_B_97_614_0;
            tmpForInput[8] = rtb_B_97_619_0;
            rtb_B_97_623_0 = -0.0;
            for (i = 0; i < 9; i++) {
                rtb_B_97_623_0 += tmpForInput[i];
            }
            B_66_38_0 = _rtP->P_445 * rtb_B_97_623_0;
            B_5_0_0_idx_3 = muDoubleScalarHypot(rtb_B_97_622_0, rtb_B_97_82_0);
            if (B_5_0_0_idx_3 > _rtP->P_446) {
                B_5_0_0_idx_3 = _rtP->P_446;
            } else {
                if (B_5_0_0_idx_3 < _rtP->P_447) {
                    B_5_0_0_idx_3 = _rtP->P_447;
                }
            }
            _rtB->B_66_40_0 = 1.0 / B_5_0_0_idx_3;
        }
        if (ssIsMajorTimeStep(S) != 0) {
            srUpdateBC(_rtDW->AutomaticGainControl_SubsysRanBC_n);
        }
    }
    rtb_B_97_623_0 = rtb_B_97_584_0 * _rtB->B_66_40_0;
    _rtB->B_97_863_0 = _rtP->P_1142 * rtb_B_97_623_0 * _rtP->P_1143;
    B_97_864_0 = _rtDW->UD_DSTATE_b;
    B_97_865_0 = _rtB->B_97_863_0 - B_97_864_0;
    B_5_0_0_idx_3 = (_rtP->P_1137 * rtb_B_97_623_0 + _rtDW->DiscreteTimeIntegrator_DSTATE_pw) + B_97_865_0;
    if (B_5_0_0_idx_3 > _rtP->P_1145) {
        _rtB->B_97_867_0 = _rtP->P_1145;
    } else if (B_5_0_0_idx_3 < _rtP->P_1146) {
        _rtB->B_97_867_0 = _rtP->P_1146;
    } else {
        _rtB->B_97_867_0 = B_5_0_0_idx_3;
    }
    B_5_0_0_idx_3 = _rtP->P_1147 * _rtB->B_97_867_0 - _rtDW->UnitDelay_DSTATE_kd;
    if (B_5_0_0_idx_3 > _rtP->P_1149) {
        B_5_0_0_idx_3 = _rtP->P_1149;
    } else {
        if (B_5_0_0_idx_3 < _rtP->P_1150) {
            B_5_0_0_idx_3 = _rtP->P_1150;
        }
    }
    _rtB->B_97_872_0 = B_5_0_0_idx_3 + _rtDW->UnitDelay_DSTATE_kd;
          {
                 _rtB->B_97_873_0 = (_rtP->P_1153[0])*_rtDW->DiscreteStateSpace_DSTATE_i[0] + (_rtP->P_1153[1])*_rtDW->DiscreteStateSpace_DSTATE_i[1];
    _rtB->B_97_873_0 += _rtP->P_1154*_rtB->B_97_872_0;
  }

   
  
    _rtB->B_97_875_0 = _rtP->P_1156 * rtb_B_97_623_0;
    B_97_876_0 = ((B_97_719_0_idx_0 * rtb_B_97_574_0 + B_97_719_0_idx_1 * rtb_B_97_805_0) + B_97_719_0_idx_2 * rtb_B_97_580_0) * 0.66666666666666663;
    tmpForInput[0] = B_97_719_0_idx_0;
    tmpForInput[1] = B_97_719_0_idx_1;
    tmpForInput[2] = B_97_719_0_idx_2;
    tmpForInput[3] = rtb_B_97_574_0;
    tmpForInput[4] = rtb_B_97_803_0;
    tmpForInput[5] = rtb_B_97_805_0;
    tmpForInput[6] = rtb_B_97_575_0;
    tmpForInput[7] = rtb_B_97_580_0;
    tmpForInput[8] = rtb_B_97_583_0;
    B_97_719_0_idx_1 = -0.0;
    for (i = 0; i < 9; i++) {
        B_97_719_0_idx_1 += tmpForInput[i];
    }
    B_97_878_0 = _rtP->P_1157 * B_97_719_0_idx_1;
    B_97_726_1_idx_0 = (_rtP->P_1158 * B_97_726_1_idx_0 >= _rtB->B_97_791_0);
    B_97_726_1_idx_1 = (_rtP->P_1158 * B_97_726_1_idx_1 >= _rtB->B_97_791_0);
    B_97_726_1_idx_2 = (_rtP->P_1158 * B_97_726_1_idx_2 >= _rtB->B_97_791_0);
    B_97_882_0[0] = !(B_97_726_1_idx_0 != 0.0);
    B_97_882_0[1] = !(B_97_726_1_idx_1 != 0.0);
    B_97_882_0[2] = !(B_97_726_1_idx_2 != 0.0);
}
_rtB->B_97_883_0 = 0.0;
_rtB->B_97_883_0 += _rtP->P_1160 * _rtX->TransferFcn4_CSTATE;
_rtB->B_97_884_0 = _rtX->Integrator_CSTATE;
_rtB->B_97_885_0 = _rtP->P_1162 * _rtB->B_97_883_0;
_rtB->B_97_886_0 = _rtX->Filter_CSTATE;
_rtB->B_97_887_0 = _rtB->B_97_885_0 - _rtB->B_97_886_0;
_rtB->B_97_888_0 = _rtP->P_1164 * _rtB->B_97_887_0;
_rtB->B_97_889_0 = 0.0;
_rtB->B_97_889_0 += _rtP->P_1166 * _rtX->TransferFcn1_CSTATE_m;
_rtB->B_97_891_0 = _rtB->B_97_889_0 - _rtB->B_97_890_0;
if (ssGetT(S) >= _rtP->P_1169) {
    _rtB->B_97_894_0 = _rtB->B_97_891_0;
} else {
    _rtB->B_97_894_0 = _rtB->B_97_893_0;
}
_rtB->B_97_895_0 = _rtX->Integrator_CSTATE_c;
_rtB->B_97_896_0 = _rtP->P_1171 * _rtB->B_97_894_0;
_rtB->B_97_897_0 = _rtX->Filter_CSTATE_e;
_rtB->B_97_898_0 = _rtB->B_97_896_0 - _rtB->B_97_897_0;
_rtB->B_97_899_0 = _rtP->P_1173 * _rtB->B_97_898_0;
_rtB->B_97_900_0 = 0.0;
_rtB->B_97_900_0 += _rtP->P_1175 * _rtX->TransferFcn2_CSTATE;
_rtB->B_97_902_0 = _rtX->Integrator_CSTATE_a;
_rtB->B_97_903_0 = _rtP->P_1178 * _rtB->B_97_900_0;
_rtB->B_97_904_0 = _rtX->Filter_CSTATE_g;
_rtB->B_97_905_0 = _rtB->B_97_903_0 - _rtB->B_97_904_0;
_rtB->B_97_906_0 = _rtP->P_1180 * _rtB->B_97_905_0;
if (ssGetT(S) >= _rtP->P_1181) {
    B_97_719_0_idx_1 = ((_rtP->P_502 * _rtB->B_97_894_0 + _rtB->B_97_895_0) + _rtB->B_97_899_0) + ((_rtP->P_1176 * _rtB->B_97_900_0 + _rtB->B_97_902_0) + _rtB->B_97_906_0);
} else {
    B_97_719_0_idx_1 = (_rtP->P_1176 * _rtB->B_97_900_0 + _rtB->B_97_902_0) + _rtB->B_97_906_0;
}
if (ssGetT(S) >= _rtP->P_1182) {
    B_97_719_0_idx_1 += (_rtP->P_501 * _rtB->B_97_883_0 + _rtB->B_97_884_0) + _rtB->B_97_888_0;
}
_rtB->B_97_914_0 = _rtP->P_1183 * B_97_719_0_idx_1;
_rtB->B_97_915_0 = 0.0;
_rtB->B_97_915_0 += _rtP->P_1185 * _rtX->TransferFcn4_CSTATE_k;
_rtB->B_97_916_0 = _rtX->Integrator_CSTATE_h;
_rtB->B_97_917_0 = _rtP->P_1187 * _rtB->B_97_915_0;
_rtB->B_97_918_0 = _rtX->Filter_CSTATE_c;
_rtB->B_97_919_0 = _rtB->B_97_917_0 - _rtB->B_97_918_0;
_rtB->B_97_920_0 = _rtP->P_1189 * _rtB->B_97_919_0;
_rtB->B_97_922_0 = _rtB->B_97_760_0 - _rtB->B_97_921_0;
if (ssGetT(S) >= _rtP->P_1192) {
    _rtB->B_97_925_0 = _rtB->B_97_922_0;
} else {
    _rtB->B_97_925_0 = _rtB->B_97_924_0;
}
_rtB->B_97_926_0 = _rtX->Integrator_CSTATE_g;
_rtB->B_97_927_0 = _rtP->P_1194 * _rtB->B_97_925_0;
_rtB->B_97_928_0 = _rtX->Filter_CSTATE_b;
_rtB->B_97_929_0 = _rtB->B_97_927_0 - _rtB->B_97_928_0;
_rtB->B_97_930_0 = _rtP->P_1196 * _rtB->B_97_929_0;
_rtB->B_97_931_0 = 0.0;
_rtB->B_97_931_0 += _rtP->P_1198 * _rtX->TransferFcn2_CSTATE_d;
_rtB->B_97_933_0 = _rtX->Integrator_CSTATE_h0;
_rtB->B_97_934_0 = _rtP->P_1201 * _rtB->B_97_931_0;
_rtB->B_97_935_0 = _rtX->Filter_CSTATE_p;
_rtB->B_97_936_0 = _rtB->B_97_934_0 - _rtB->B_97_935_0;
_rtB->B_97_937_0 = _rtP->P_1203 * _rtB->B_97_936_0;
if (ssGetT(S) >= _rtP->P_1204) {
    B_97_719_0_idx_1 = ((_rtP->P_536 * _rtB->B_97_925_0 + _rtB->B_97_926_0) + _rtB->B_97_930_0) + ((_rtP->P_1199 * _rtB->B_97_931_0 + _rtB->B_97_933_0) + _rtB->B_97_937_0);
} else {
    B_97_719_0_idx_1 = (_rtP->P_1199 * _rtB->B_97_931_0 + _rtB->B_97_933_0) + _rtB->B_97_937_0;
}
if (ssGetT(S) >= _rtP->P_1205) {
    B_97_719_0_idx_1 += (_rtP->P_535 * _rtB->B_97_915_0 + _rtB->B_97_916_0) + _rtB->B_97_920_0;
}
_rtB->B_97_945_0 = _rtP->P_1206 * B_97_719_0_idx_1;
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 946, SS_CALL_MDL_OUTPUTS);
   
  
_rtB->B_97_949_0 = _rtP->P_1207 * _rtB->B_97_900_0;
_rtB->B_97_950_0 = _rtP->P_1208 * _rtB->B_97_894_0;
_rtB->B_97_951_0 = _rtP->P_1209 * _rtB->B_97_883_0;
if (ssGetT(S) >= _rtP->P_1211) {
    _rtB->B_97_953_0 = _rtB->B_97_891_0;
} else {
    _rtB->B_97_953_0 = _rtB->B_97_952_0;
}
if (ssGetT(S) >= _rtP->P_1213) {
    _rtB->B_97_955_0 = _rtB->B_97_891_0;
} else {
    _rtB->B_97_955_0 = _rtB->B_97_954_0;
}
_rtB->B_97_957_0 = _rtB->B_97_956_0 + _rtB->B_97_914_0;
if (ssIsMajorTimeStep(S) != 0) {
    _rtDW->Saturation_MODE = _rtB->B_97_957_0 >= _rtP->P_1215 ? 1 : _rtB->B_97_957_0 > _rtP->P_1216 ? 0 : -1;
}
B_97_719_0_idx_1 = _rtDW->Saturation_MODE == 1 ? _rtP->P_1215 : _rtDW->Saturation_MODE == -1 ? _rtP->P_1216 : _rtB->B_97_957_0;
_rtB->B_97_960_0 = _rtB->B_97_959_0 + _rtB->B_97_945_0;
if (ssIsMajorTimeStep(S) != 0) {
    _rtDW->Saturation1_MODE = _rtB->B_97_960_0 >= _rtP->P_1218 ? 1 : _rtB->B_97_960_0 > _rtP->P_1219 ? 0 : -1;
}
rtb_B_97_623_0 = _rtDW->Saturation1_MODE == 1 ? _rtP->P_1218 : _rtDW->Saturation1_MODE == -1 ? _rtP->P_1219 : _rtB->B_97_960_0;
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_964_0 = muDoubleScalarMod(_rtDW->DiscreteTimeIntegrator_DSTATE_o, _rtB->B_97_963_0);
    B_97_965_0_idx_0 = _rtP->P_1223 * B_97_291_1_idx_2;
    B_97_965_0_idx_1 = _rtP->P_1223 * B_97_328_1_idx_1;
    B_97_965_0_idx_2 = _rtP->P_1223 * B_97_328_1_idx_2;
    rtb_B_97_966_0 = muDoubleScalarSin(_rtB->B_97_964_0);
    rtb_B_97_967_0 = muDoubleScalarCos(_rtB->B_97_964_0);
    rtb_B_97_972_0 = (0.0 - rtb_B_97_966_0 * _rtB->B_97_968_0) - rtb_B_97_967_0 * _rtB->B_97_970_0;
    rtb_B_97_975_0 = rtb_B_97_966_0 * _rtB->B_97_970_0 - rtb_B_97_967_0 * _rtB->B_97_968_0;
    rtb_B_97_976_0 = (0.0 - rtb_B_97_972_0) - rtb_B_97_966_0;
    rtb_B_97_977_0 = (0.0 - rtb_B_97_975_0) - rtb_B_97_967_0;
    _rtB->B_97_978_0 = ((B_97_965_0_idx_0 * rtb_B_97_966_0 + B_97_965_0_idx_1 * rtb_B_97_972_0) + B_97_965_0_idx_2 * rtb_B_97_976_0) * 0.66666666666666663;
    _rtB->B_97_979_0 = ((B_97_965_0_idx_0 * rtb_B_97_967_0 + B_97_965_0_idx_1 * rtb_B_97_975_0) + B_97_965_0_idx_2 * rtb_B_97_977_0) * 0.66666666666666663;
    B_97_726_1_idx_0 = _rtP->P_1226 * B_97_321_0_idx_0 * _rtP->P_1227;
    B_97_726_1_idx_1 = _rtP->P_1226 * B_97_321_0_idx_1 * _rtP->P_1227;
    B_97_726_1_idx_2 = _rtP->P_1226 * B_97_321_0_idx_2 * _rtP->P_1227;
    rtb_B_97_986_0 = (0.0 - rtb_B_97_966_0 * _rtB->B_97_982_0) - rtb_B_97_967_0 * _rtB->B_97_984_0;
    rtb_B_97_989_0 = rtb_B_97_966_0 * _rtB->B_97_984_0 - rtb_B_97_967_0 * _rtB->B_97_982_0;
    rtb_B_97_990_0 = (0.0 - rtb_B_97_986_0) - rtb_B_97_966_0;
    rtb_B_97_991_0 = (0.0 - rtb_B_97_989_0) - rtb_B_97_967_0;
    _rtB->B_97_992_0 = ((B_97_726_1_idx_0 * rtb_B_97_966_0 + B_97_726_1_idx_1 * rtb_B_97_986_0) + B_97_726_1_idx_2 * rtb_B_97_990_0) * 0.66666666666666663;
    _rtB->B_97_993_0 = ((B_97_726_1_idx_0 * rtb_B_97_967_0 + B_97_726_1_idx_1 * rtb_B_97_989_0) + B_97_726_1_idx_2 * rtb_B_97_991_0) * 0.66666666666666663;
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 994, SS_CALL_MDL_OUTPUTS);
   
  
}
_rtB->B_97_996_0 = (_rtB->B_97_657_0[0] - B_97_719_0_idx_1) * _rtP->P_1230;
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_998_0 = _rtP->P_1232 * _rtB->B_97_996_0;
    _rtB->B_97_999_0 = _rtP->P_1233 * _rtB->B_97_998_0 + _rtDW->Integrator_DSTATE;
    B_5_0_0_idx_3 = _rtP->P_1231 * _rtB->B_97_996_0 + _rtB->B_97_999_0;
    if (B_5_0_0_idx_3 > _rtP->P_1235) {
        _rtB->B_97_1001_0 = _rtP->P_1235;
    } else if (B_5_0_0_idx_3 < _rtP->P_1236) {
        _rtB->B_97_1001_0 = _rtP->P_1236;
    } else {
        _rtB->B_97_1001_0 = B_5_0_0_idx_3;
    }
}
_rtB->B_97_1003_0 = (_rtB->B_97_657_0[1] - rtb_B_97_623_0) * _rtP->P_1237;
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_1005_0 = _rtP->P_1239 * _rtB->B_97_1003_0;
    _rtB->B_97_1006_0 = _rtP->P_1240 * _rtB->B_97_1005_0 + _rtDW->Integrator_DSTATE_l;
    B_5_0_0_idx_3 = _rtP->P_1238 * _rtB->B_97_1003_0 + _rtB->B_97_1006_0;
    if (B_5_0_0_idx_3 > _rtP->P_1242) {
        _rtB->B_97_1008_0 = _rtP->P_1242;
    } else if (B_5_0_0_idx_3 < _rtP->P_1243) {
        _rtB->B_97_1008_0 = _rtP->P_1243;
    } else {
        _rtB->B_97_1008_0 = B_5_0_0_idx_3;
    }
    _rtB->B_97_1009_0[0] = _rtB->B_97_1001_0 - _rtB->B_97_992_0;
    _rtB->B_97_1009_0[1] = _rtB->B_97_1008_0 - _rtB->B_97_993_0;
    B_5_0_0_idx_3 = _rtP->P_1244 * _rtB->B_97_1009_0[0] + _rtDW->Integrator_DSTATE_n[0];
    if (B_5_0_0_idx_3 > _rtP->P_1247) {
        _rtB->B_97_1013_0[0] = _rtP->P_1247;
    } else if (B_5_0_0_idx_3 < _rtP->P_1248) {
        _rtB->B_97_1013_0[0] = _rtP->P_1248;
    } else {
        _rtB->B_97_1013_0[0] = B_5_0_0_idx_3;
    }
    B_5_0_0_idx_3 = _rtP->P_1244 * _rtB->B_97_1009_0[1] + _rtDW->Integrator_DSTATE_n[1];
    if (B_5_0_0_idx_3 > _rtP->P_1247) {
        _rtB->B_97_1013_0[1] = _rtP->P_1247;
    } else if (B_5_0_0_idx_3 < _rtP->P_1248) {
        _rtB->B_97_1013_0[1] = _rtP->P_1248;
    } else {
        _rtB->B_97_1013_0[1] = B_5_0_0_idx_3;
    }
    _rtB->B_97_1016_0 = (_rtP->P_1249 * _rtB->B_97_1001_0 + _rtB->B_97_978_0) - _rtP->P_1250 * _rtB->B_97_1008_0;
    _rtB->B_97_1019_0 = (_rtP->P_1251 * _rtB->B_97_1008_0 + _rtB->B_97_979_0) + _rtP->P_1252 * _rtB->B_97_1001_0;
    _rtB->B_97_1021_0[0] = (_rtB->B_97_1013_0[0] + _rtB->B_97_1016_0) * _rtP->P_1253;
    _rtB->B_97_1021_0[1] = (_rtB->B_97_1013_0[1] + _rtB->B_97_1019_0) * _rtP->P_1253;
    if (_rtB->B_97_1021_0[0] > _rtP->P_1254) {
        _rtB->B_97_1022_0[0] = _rtP->P_1254;
    } else if (_rtB->B_97_1021_0[0] < _rtP->P_1255) {
        _rtB->B_97_1022_0[0] = _rtP->P_1255;
    } else {
        _rtB->B_97_1022_0[0] = _rtB->B_97_1021_0[0];
    }
    if (_rtB->B_97_1021_0[1] > _rtP->P_1254) {
        _rtB->B_97_1022_0[1] = _rtP->P_1254;
    } else if (_rtB->B_97_1021_0[1] < _rtP->P_1255) {
        _rtB->B_97_1022_0[1] = _rtP->P_1255;
    } else {
        _rtB->B_97_1022_0[1] = _rtB->B_97_1021_0[1];
    }
    B_5_0_0_idx_3 = _rtB->B_97_1022_0[0] / _rtB->B_97_1026_0;
    rtb_B_97_176_0 = _rtB->B_97_1022_0[1] / _rtB->B_97_1026_0;
    B_97_719_0_idx_1 = muDoubleScalarAtan2(rtb_B_97_176_0, B_5_0_0_idx_3);
    rtb_B_97_623_0 = muDoubleScalarHypot(B_5_0_0_idx_3, rtb_B_97_176_0);
    if (rtb_B_97_623_0 > _rtP->P_1259) {
        rtb_B_97_623_0 = _rtP->P_1259;
    } else {
        if (rtb_B_97_623_0 < _rtP->P_1260) {
            rtb_B_97_623_0 = _rtP->P_1260;
        }
    }
    _rtB->B_97_1037_0[0] = muDoubleScalarSin((((_rtB->B_97_964_0 + _rtB->B_97_1031_0[0]) + _rtB->B_97_1032_0) + _rtB->B_97_1033_0) + B_97_719_0_idx_1) * rtb_B_97_623_0;
    _rtB->B_97_1037_0[1] = muDoubleScalarSin((((_rtB->B_97_964_0 + _rtB->B_97_1031_0[1]) + _rtB->B_97_1032_0) + _rtB->B_97_1033_0) + B_97_719_0_idx_1) * rtb_B_97_623_0;
    _rtB->B_97_1037_0[2] = muDoubleScalarSin((((_rtB->B_97_964_0 + _rtB->B_97_1031_0[2]) + _rtB->B_97_1032_0) + _rtB->B_97_1033_0) + B_97_719_0_idx_1) * rtb_B_97_623_0;
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 1038, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 1039, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 1040, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 1041, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 1042, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 1043, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 1044, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 1045, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 1046, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 1047, SS_CALL_MDL_OUTPUTS);
   
  
    _rtB->B_97_1048_0[0] = _rtP->P_1264 * _rtB->B_97_1009_0[0];
    _rtB->B_97_1048_0[1] = _rtP->P_1264 * _rtB->B_97_1009_0[1];
}
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 1049, SS_CALL_MDL_OUTPUTS);
   
  
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    B_97_719_0_idx_1 = (0.0 - rtb_B_97_966_0 * _rtB->B_97_1050_0) - rtb_B_97_967_0 * _rtB->B_97_1052_0;
    rtb_B_97_623_0 = rtb_B_97_966_0 * _rtB->B_97_1052_0 - rtb_B_97_967_0 * _rtB->B_97_1050_0;
    rtb_B_97_574_0 = (0.0 - B_97_719_0_idx_1) - rtb_B_97_966_0;
    rtb_B_97_803_0 = (0.0 - rtb_B_97_623_0) - rtb_B_97_967_0;
    _rtB->B_97_1060_0 = ((B_97_965_0_idx_0 * rtb_B_97_967_0 + B_97_965_0_idx_1 * rtb_B_97_623_0) + B_97_965_0_idx_2 * rtb_B_97_803_0) * 0.66666666666666663;
    if (_rtDW->Integ4_SYSTEM_ENABLE_hq != 0) {
        _rtB->B_97_1061_0 = _rtDW->Integ4_DSTATE_bi;
    } else {
        _rtB->B_97_1061_0 = _rtP->P_1267 * _rtB->B_97_1060_0 + _rtDW->Integ4_DSTATE_bi;
    }
    rtb_B_97_805_0 = _rtDW->UnitDelay_DSTATE_dd;
    if (_rtDW->UnitDelay_DSTATE_dd > _rtP->P_1270) {
        B_5_0_0_idx_3 = _rtP->P_1270;
    } else if (_rtDW->UnitDelay_DSTATE_dd < _rtP->P_1271) {
        B_5_0_0_idx_3 = _rtP->P_1271;
    } else {
        B_5_0_0_idx_3 = _rtDW->UnitDelay_DSTATE_dd;
    }
    rtb_B_97_575_0 = 1.0 / B_5_0_0_idx_3 / 1.0e-5;
    rtb_B_97_580_0 = muDoubleScalarCeil(rtb_B_97_575_0);
    _rtB->B_97_1066_0 = _rtP->P_1272 * rtb_B_97_580_0;
                
          /* Level2 S-Function Block: '<S1155>/B_82_190' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 1067, SS_CALL_MDL_OUTPUTS);

   
  
    rtb_B_97_82_0 = ssGetTaskTime(S, 2);
    if (rtb_B_97_82_0 < _rtP->P_1274) {
        rtb_B_97_82_0 = _rtP->P_1275;
    } else {
        rtb_B_97_82_0 = _rtP->P_1276;
    }
    if (rtb_B_97_82_0 >= _rtP->P_1278) {
        rtb_B_97_82_0 = rtb_B_97_575_0 - rtb_B_97_580_0;
        rtb_B_97_575_0 = ((_rtB->B_97_1060_0 - _rtDW->UnitDelay_DSTATE_mf) * rtb_B_97_82_0 * _rtP->P_534 + _rtB->B_97_1060_0) * (rtb_B_97_82_0 / rtb_B_97_575_0) + (_rtB->B_97_1061_0 - _rtB->B_97_1067_0) * rtb_B_97_805_0;
    } else {
        rtb_B_97_575_0 = _rtP->P_1277;
    }
    _rtB->B_97_1074_0 = _rtB->B_97_1073_0;
    i = ssIsSampleHit(S, 2, 0);
    if ((i != 0) && (ssIsMajorTimeStep(S) != 0)) {
        if (_rtB->B_97_1074_0 > 0.0) {
            if (!_rtDW->AutomaticGainControl_MODE) {
                if (ssGetTaskTime(S, 2) != ssGetTStart(S)) {
                    ssSetBlockStateForSolverChangedAtMajorStep(S);
                }
  
                _rtDW->Integ4_SYSTEM_ENABLE_f0v = 1U;
                _rtDW->Integ4_SYSTEM_ENABLE_ki = 1U;
                _rtDW->AutomaticGainControl_MODE = true;
            }
        } else {
            if (_rtDW->AutomaticGainControl_MODE) {
                ssSetBlockStateForSolverChangedAtMajorStep(S);
  
                _rtDW->AutomaticGainControl_MODE = false;
            }
        }
    }
    if (_rtDW->AutomaticGainControl_MODE) {
        i = ssIsSampleHit(S, 2, 0);
        if (i != 0) {
            rtb_B_97_580_0 = (0.0 - rtb_B_97_966_0 * _rtB->B_93_0_0) - rtb_B_97_967_0 * _rtB->B_93_2_0;
            rtb_B_97_583_0 = rtb_B_97_966_0 * _rtB->B_93_2_0 - rtb_B_97_967_0 * _rtB->B_93_0_0;
            rtb_B_97_584_0 = (0.0 - rtb_B_97_580_0) - rtb_B_97_966_0;
            rtb_B_97_585_0 = (0.0 - rtb_B_97_583_0) - rtb_B_97_967_0;
            _rtB->B_93_10_0 = ((B_97_965_0_idx_0 * rtb_B_97_966_0 + B_97_965_0_idx_1 * rtb_B_97_580_0) + B_97_965_0_idx_2 * rtb_B_97_584_0) * 0.66666666666666663;
            if (_rtDW->Integ4_SYSTEM_ENABLE_f0v != 0) {
                _rtB->B_93_11_0 = _rtDW->Integ4_DSTATE_pa;
            } else {
                _rtB->B_93_11_0 = _rtP->P_508 * _rtB->B_93_10_0 + _rtDW->Integ4_DSTATE_pa;
            }
            if (rtb_B_97_805_0 > _rtP->P_510) {
                B_97_321_0_idx_0 = _rtP->P_510;
            } else if (rtb_B_97_805_0 < _rtP->P_511) {
                B_97_321_0_idx_0 = _rtP->P_511;
            } else {
                B_97_321_0_idx_0 = rtb_B_97_805_0;
            }
            rtb_B_97_613_0 = 1.0 / B_97_321_0_idx_0 / 1.0e-5;
            rtb_B_97_614_0 = muDoubleScalarCeil(rtb_B_97_613_0);
            _rtB->B_93_15_0 = _rtP->P_512 * rtb_B_97_614_0;
                
          /* Level2 S-Function Block: '<S1150>/B_78_0' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 93, 16, SS_CALL_MDL_OUTPUTS);

   
  
            rtb_B_97_82_0 = ssGetTaskTime(S, 2);
            if (rtb_B_97_82_0 < _rtP->P_514) {
                rtb_B_97_82_0 = _rtP->P_515;
            } else {
                rtb_B_97_82_0 = _rtP->P_516;
            }
            if (rtb_B_97_82_0 >= _rtP->P_518) {
                rtb_B_97_82_0 = rtb_B_97_613_0 - rtb_B_97_614_0;
                rtb_B_97_613_0 = ((_rtB->B_93_10_0 - _rtDW->UnitDelay_DSTATE_eer) * rtb_B_97_82_0 * _rtP->P_503 + _rtB->B_93_10_0) * (rtb_B_97_82_0 / rtb_B_97_613_0) + (_rtB->B_93_11_0 - _rtB->B_93_16_0) * rtb_B_97_805_0;
            } else {
                rtb_B_97_613_0 = _rtP->P_517;
            }
            _rtB->B_93_22_0 = ((B_97_965_0_idx_0 * rtb_B_97_967_0 + B_97_965_0_idx_1 * rtb_B_97_583_0) + B_97_965_0_idx_2 * rtb_B_97_585_0) * 0.66666666666666663;
            if (_rtDW->Integ4_SYSTEM_ENABLE_ki != 0) {
                _rtB->B_93_23_0 = _rtDW->Integ4_DSTATE_jc;
            } else {
                _rtB->B_93_23_0 = _rtP->P_519 * _rtB->B_93_22_0 + _rtDW->Integ4_DSTATE_jc;
            }
            if (rtb_B_97_805_0 > _rtP->P_521) {
                B_97_321_0_idx_0 = _rtP->P_521;
            } else if (rtb_B_97_805_0 < _rtP->P_522) {
                B_97_321_0_idx_0 = _rtP->P_522;
            } else {
                B_97_321_0_idx_0 = rtb_B_97_805_0;
            }
            rtb_B_97_614_0 = 1.0 / B_97_321_0_idx_0 / 1.0e-5;
            rtb_B_97_619_0 = muDoubleScalarCeil(rtb_B_97_614_0);
            _rtB->B_93_27_0 = _rtP->P_523 * rtb_B_97_619_0;
                
          /* Level2 S-Function Block: '<S1152>/B_78_2' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 93, 28, SS_CALL_MDL_OUTPUTS);

   
  
            rtb_B_97_82_0 = ssGetTaskTime(S, 2);
            if (rtb_B_97_82_0 < _rtP->P_525) {
                rtb_B_97_82_0 = _rtP->P_526;
            } else {
                rtb_B_97_82_0 = _rtP->P_527;
            }
            if (rtb_B_97_82_0 >= _rtP->P_529) {
                rtb_B_97_82_0 = rtb_B_97_614_0 - rtb_B_97_619_0;
                rtb_B_97_82_0 = ((_rtB->B_93_22_0 - _rtDW->UnitDelay_DSTATE_fv) * rtb_B_97_82_0 * _rtP->P_504 + _rtB->B_93_22_0) * (rtb_B_97_82_0 / rtb_B_97_614_0) + (_rtB->B_93_23_0 - _rtB->B_93_28_0) * rtb_B_97_805_0;
            } else {
                rtb_B_97_82_0 = _rtP->P_528;
            }
            B_93_36_0 = _rtP->P_530 * muDoubleScalarAtan2(rtb_B_97_82_0, rtb_B_97_613_0);
            tmpForInput[0] = B_97_965_0_idx_0;
            tmpForInput[1] = B_97_965_0_idx_1;
            tmpForInput[2] = B_97_965_0_idx_2;
            tmpForInput[3] = rtb_B_97_966_0;
            tmpForInput[4] = rtb_B_97_967_0;
            tmpForInput[5] = rtb_B_97_580_0;
            tmpForInput[6] = rtb_B_97_583_0;
            tmpForInput[7] = rtb_B_97_584_0;
            tmpForInput[8] = rtb_B_97_585_0;
            B_97_321_0_idx_0 = -0.0;
            for (i = 0; i < 9; i++) {
                B_97_321_0_idx_0 += tmpForInput[i];
            }
            B_93_38_0 = _rtP->P_531 * B_97_321_0_idx_0;
            B_5_0_0_idx_3 = muDoubleScalarHypot(rtb_B_97_613_0, rtb_B_97_82_0);
            if (B_5_0_0_idx_3 > _rtP->P_532) {
                B_5_0_0_idx_3 = _rtP->P_532;
            } else {
                if (B_5_0_0_idx_3 < _rtP->P_533) {
                    B_5_0_0_idx_3 = _rtP->P_533;
                }
            }
            _rtB->B_93_40_0 = 1.0 / B_5_0_0_idx_3;
        }
        if (ssIsMajorTimeStep(S) != 0) {
            srUpdateBC(_rtDW->AutomaticGainControl_SubsysRanBC);
        }
    }
    B_97_321_0_idx_0 = rtb_B_97_575_0 * _rtB->B_93_40_0;
    _rtB->B_97_1081_0 = _rtP->P_1285 * B_97_321_0_idx_0 * _rtP->P_1286;
    B_97_1082_0 = _rtDW->UD_DSTATE_a;
    B_97_1083_0 = _rtB->B_97_1081_0 - B_97_1082_0;
    B_5_0_0_idx_3 = (_rtP->P_1280 * B_97_321_0_idx_0 + _rtDW->DiscreteTimeIntegrator_DSTATE_jp) + B_97_1083_0;
    if (B_5_0_0_idx_3 > _rtP->P_1288) {
        _rtB->B_97_1085_0 = _rtP->P_1288;
    } else if (B_5_0_0_idx_3 < _rtP->P_1289) {
        _rtB->B_97_1085_0 = _rtP->P_1289;
    } else {
        _rtB->B_97_1085_0 = B_5_0_0_idx_3;
    }
    B_5_0_0_idx_3 = _rtP->P_1290 * _rtB->B_97_1085_0 - _rtDW->UnitDelay_DSTATE_kv;
    if (B_5_0_0_idx_3 > _rtP->P_1292) {
        B_5_0_0_idx_3 = _rtP->P_1292;
    } else {
        if (B_5_0_0_idx_3 < _rtP->P_1293) {
            B_5_0_0_idx_3 = _rtP->P_1293;
        }
    }
    _rtB->B_97_1090_0 = B_5_0_0_idx_3 + _rtDW->UnitDelay_DSTATE_kv;
          {
                 _rtB->B_97_1091_0 = (_rtP->P_1296[0])*_rtDW->DiscreteStateSpace_DSTATE_p[0] + (_rtP->P_1296[1])*_rtDW->DiscreteStateSpace_DSTATE_p[1];
    _rtB->B_97_1091_0 += _rtP->P_1297*_rtB->B_97_1090_0;
  }

   
  
    _rtB->B_97_1093_0 = _rtP->P_1299 * B_97_321_0_idx_0;
    B_97_1094_0 = ((B_97_965_0_idx_0 * rtb_B_97_966_0 + B_97_965_0_idx_1 * B_97_719_0_idx_1) + B_97_965_0_idx_2 * rtb_B_97_574_0) * 0.66666666666666663;
    tmpForInput[0] = B_97_965_0_idx_0;
    tmpForInput[1] = B_97_965_0_idx_1;
    tmpForInput[2] = B_97_965_0_idx_2;
    tmpForInput[3] = rtb_B_97_966_0;
    tmpForInput[4] = rtb_B_97_967_0;
    tmpForInput[5] = B_97_719_0_idx_1;
    tmpForInput[6] = rtb_B_97_623_0;
    tmpForInput[7] = rtb_B_97_574_0;
    tmpForInput[8] = rtb_B_97_803_0;
    B_97_321_0_idx_0 = -0.0;
    for (i = 0; i < 9; i++) {
        B_97_321_0_idx_0 += tmpForInput[i];
    }
    B_97_1096_0 = _rtP->P_1300 * B_97_321_0_idx_0;
    tmpForInput[0] = B_97_965_0_idx_0;
    tmpForInput[1] = B_97_965_0_idx_1;
    tmpForInput[2] = B_97_965_0_idx_2;
    tmpForInput[3] = rtb_B_97_966_0;
    tmpForInput[4] = rtb_B_97_967_0;
    tmpForInput[5] = rtb_B_97_972_0;
    tmpForInput[6] = rtb_B_97_975_0;
    tmpForInput[7] = rtb_B_97_976_0;
    tmpForInput[8] = rtb_B_97_977_0;
    rtb_B_97_972_0 = -0.0;
    for (i = 0; i < 9; i++) {
        rtb_B_97_972_0 += tmpForInput[i];
    }
    B_97_1098_0 = _rtP->P_1301 * rtb_B_97_972_0;
    tmpForInput[0] = B_97_726_1_idx_0;
    tmpForInput[1] = B_97_726_1_idx_1;
    tmpForInput[2] = B_97_726_1_idx_2;
    tmpForInput[3] = rtb_B_97_966_0;
    tmpForInput[4] = rtb_B_97_967_0;
    tmpForInput[5] = rtb_B_97_986_0;
    tmpForInput[6] = rtb_B_97_989_0;
    tmpForInput[7] = rtb_B_97_990_0;
    tmpForInput[8] = rtb_B_97_991_0;
    rtb_B_97_966_0 = -0.0;
    for (i = 0; i < 9; i++) {
        rtb_B_97_966_0 += tmpForInput[i];
    }
    B_97_1100_0 = _rtP->P_1302 * rtb_B_97_966_0;
    rtb_B_97_966_0 = look1_pbinlxpw(muDoubleScalarRem(ssGetTaskTime(S, 2) + _rtP->P_1303, _rtP->P_1304) * _rtP->P_1305, _rtP->P_1307, _rtP->P_1306, &_rtDW->m_bpIndex, 2U);
    B_97_1109_0 = rtb_B_97_966_0 - _rtP->P_1308;
}
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    _rtB->B_97_1111_0 = _rtP->P_1309;
}
rtb_B_97_966_0 = ssGetT(S) + _rtB->B_97_1111_0;
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    _rtB->B_97_1113_0 = _rtP->P_1310;
}
rtb_B_97_966_0 = look1_pbinlxpw(_rtP->P_1311 * muDoubleScalarRem(rtb_B_97_966_0, _rtB->B_97_1113_0), _rtP->P_1313, _rtP->P_1312, &_rtDW->m_bpIndex_g, 2U);
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    _rtB->B_97_1117_0 = _rtP->P_1314;
}
B_97_1118_0 = rtb_B_97_966_0 - _rtB->B_97_1117_0;
i = ssIsSampleHit(S, 4, 0);
if (i != 0) {
    _rtB->B_97_1119_0 = _rtDW->UnitDelay_DSTATE_ef;
}
_rtB->B_97_1122_0 = _rtP->P_1316 * _rtB->B_97_931_0;
_rtB->B_97_1123_0 = _rtP->P_1317 * _rtB->B_97_925_0;
_rtB->B_97_1124_0 = _rtP->P_1318 * _rtB->B_97_915_0;
if (ssGetT(S) >= _rtP->P_1320) {
    _rtB->B_97_1126_0 = _rtB->B_97_922_0;
} else {
    _rtB->B_97_1126_0 = _rtB->B_97_1125_0;
}
if (ssGetT(S) >= _rtP->P_1322) {
    _rtB->B_97_1128_0 = _rtB->B_97_922_0;
} else {
    _rtB->B_97_1128_0 = _rtB->B_97_1127_0;
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    if (_rtB->B_97_1169_0 >= _rtP->P_1329) {
        rtb_B_97_176_0 = _rtB->Saturation_i.B_3_10_0;
    } else {
        rtb_B_97_176_0 = _rtB->B_97_1170_0[1];
    }
    B_97_726_1_idx_0 = _rtP->P_1330[0] * _rtB->B_97_24_0[1] * rtb_B_97_176_0;
    B_97_726_1_idx_1 = _rtP->P_1330[1] * _rtB->B_97_24_0[2] * rtb_B_97_176_0;
    B_97_726_1_idx_2 = _rtP->P_1330[2] * _rtB->B_97_24_0[3] * rtb_B_97_176_0;
    i = ssIsSampleHit(S, 2, 0);
    if (i != 0) {
        rtb_B_97_1177_0[0] = _rtB->B_97_116_0[25];
        rtb_B_97_1177_0[1] = _rtB->B_97_116_0[26];
    }
    rtb_B_97_1177_0[0] *= _rtP->P_1331;
    rtb_B_97_1177_0[1] *= _rtP->P_1331;
    _rtB->B_97_1179_0 = ((2.0 * rtb_B_97_1177_0[0] + rtb_B_97_1177_0[1]) * rtb_B_97_5_1 + 1.7320508075688772 * rtb_B_97_1177_0[1] * rtb_B_97_5_0) * 0.33333333333333331;
    _rtB->B_97_1180_0 = ((2.0 * rtb_B_97_1177_0[0] + rtb_B_97_1177_0[1]) * rtb_B_97_5_0 + -1.7320508075688772 * rtb_B_97_1177_0[1] * rtb_B_97_5_1) * 0.33333333333333331;
    rtb_B_97_1194_0[11] = _rtP->P_1338[11] * _rtB->B_97_1179_0;
    rtb_B_97_1194_0[12] = _rtP->P_1338[12] * _rtB->B_97_1180_0;
    rtb_B_97_67_0 = muDoubleScalarPower(rtb_B_97_1194_0[12], 2.0) + muDoubleScalarPower(rtb_B_97_1194_0[11], 2.0);
    if (rtb_B_97_67_0 < 0.0) {
        _rtB->B_97_1196_0 = -muDoubleScalarSqrt(-rtb_B_97_67_0);
    } else {
        _rtB->B_97_1196_0 = muDoubleScalarSqrt(rtb_B_97_67_0);
    }
          {
                 B_97_1197_0 = _rtP->P_1341*_rtDW->DiscreteStateSpace_DSTATE_n;
    B_97_1197_0 += _rtP->P_1342*_rtB->B_97_1196_0;
  }

   
  
    _rtB->B_97_1199_0 = ((_rtB->B_97_1163_0 + _rtB->B_97_1164_0) - B_97_1197_0) - _rtDW->UnitDelay1_DSTATE_j5g;
          {
                 B_97_1200_0 = _rtP->P_1348*_rtB->B_97_1199_0;
  }

   
  
    _rtB->B_97_1201_0 = _rtP->P_1350 * B_97_1200_0;
          {
                 B_97_1202_0 = _rtP->P_1353*_rtDW->DiscreteStateSpace_DSTATE_m;
    B_97_1202_0 += _rtP->P_1354*_rtB->B_97_1201_0;
  }

   
  
    rtb_B_97_1173_0 = (B_97_1202_0 > _rtB->B_97_1203_0);
    rtb_B_97_5_0 = _rtP->P_1357 * B_97_1197_0 + _rtB->B_97_1210_0;
    rtb_B_97_5_1 = (B_97_1202_0 < rtb_B_97_5_0);
    _rtB->B_97_1223_0 = (((real_T)((rtb_B_97_1173_0 != 0.0) && (rtb_B_97_5_1 != 0.0)) * B_97_1202_0 + (real_T)!(rtb_B_97_1173_0 != 0.0) * _rtB->B_97_1203_0) + (real_T)!(rtb_B_97_5_1 != 0.0) * rtb_B_97_5_0) * _rtP->P_1361;
          {
                 _rtB->B_97_1224_0 = _rtP->P_1364*_rtDW->DiscreteStateSpace_DSTATE_k;
    _rtB->B_97_1224_0 += _rtP->P_1365*_rtB->B_97_1223_0;
  }

   
  
    _rtB->B_97_1225_0 = _rtP->P_1367 * _rtB->B_97_1224_0;
          {
                 _rtB->B_97_1226_0 = (_rtP->P_1370)*_rtDW->DiscreteStateSpace_DSTATE_h;
    _rtB->B_97_1226_0 += _rtP->P_1371*_rtB->B_97_1225_0;
  }

   
  
    rtb_B_97_1173_0 = _rtDW->UnitDelay2_DSTATE_n;
    rtb_B_97_67_0 = muDoubleScalarPower(rtb_B_97_1194_0[12], 2.0) + muDoubleScalarPower(rtb_B_97_1194_0[11], 2.0);
    if (rtb_B_97_67_0 < 0.0) {
        B_97_1228_0 = -muDoubleScalarSqrt(-rtb_B_97_67_0);
    } else {
        B_97_1228_0 = muDoubleScalarSqrt(rtb_B_97_67_0);
    }
}
_rtB->B_97_1229_0 = 0.0;
_rtB->B_97_1229_0 += _rtP->P_1375[0] * _rtX->CONTROLSYSTEM_CSTATE[0];
_rtB->B_97_1229_0 += _rtP->P_1375[1] * _rtX->CONTROLSYSTEM_CSTATE[1];
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->ENGINETd_PWORK.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->ENGINETd_PWORK.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_1376;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        B_97_1230_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->ENGINETd_IWORK.CircularBufSize,
        &_rtDW->ENGINETd_IWORK.Last,
        _rtDW->ENGINETd_IWORK.Tail,
        _rtDW->ENGINETd_IWORK.Head,
        _rtP->P_1377,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
_rtB->B_97_1231_0 = _rtP->P_1378 * _rtB->B_97_1229_0;
/* Limited  Integrator  */
if (ssIsMajorTimeStep(S) != 0) {
    if (_rtX->Integrator_CSTATE_g1 >= _rtP->P_1380) {
        _rtX->Integrator_CSTATE_g1 = _rtP->P_1380;
    } else {
        if (_rtX->Integrator_CSTATE_g1 <= _rtP->P_1381) {
            _rtX->Integrator_CSTATE_g1 = _rtP->P_1381;
        }
    }
}
_rtB->B_97_1232_0 = _rtX->Integrator_CSTATE_g1;
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    rtb_B_97_1177_0[0] = _rtP->P_1383[0] * rtb_B_97_25_0[0] * _rtB->B_97_24_0[1];
    rtb_B_97_1177_0[1] = _rtP->P_1383[1] * rtb_B_97_25_0[1] * _rtB->B_97_24_0[0];
    rtb_B_97_1236_0 = rtb_B_97_1177_0[0] + rtb_B_97_1177_0[1];
    i = ssIsSampleHit(S, 2, 0);
    if (i != 0) {
        _rtB->B_97_1243_0[0] = rtb_B_97_1_0;
        _rtB->B_97_1243_0[1] = _rtP->P_1382 * rtb_B_97_17_0;
        _rtB->B_97_1243_0[2] = rtb_B_97_17_0 * rtb_B_97_1236_0 * _rtP->P_1384;
        _rtB->B_97_1243_0[3] = rtb_B_97_16_0;
        _rtB->B_97_1243_0[4] = _rtP->P_1387 * muDoubleScalarRem(_rtDW->DiscreteTimeIntegrator2_DSTATE, 6.2831853071795862);
        _rtB->B_97_1243_0[5] = _rtP->P_1388 * rtb_B_97_1236_0;
    }
}
_rtB->B_97_1244_0 = B_97_1230_0 * _rtB->B_97_1243_0[1];
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_1246_0 = _rtB->B_97_1245_0 - _rtB->B_97_1243_0[1];
}
_rtB->B_97_1247_0 = 0.0;
_rtB->B_97_1247_0 += _rtP->P_1391 * _rtX->TF1_CSTATE;
_rtB->B_97_1247_0 += _rtP->P_1392 * _rtB->B_97_1231_0;
_rtB->B_97_1248_0 = 0.0;
_rtB->B_97_1248_0 += _rtP->P_1394 * _rtX->TF2_CSTATE;
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    B_97_1249_0 = muDoubleScalarHypot(_rtB->B_97_1179_0, _rtB->B_97_1180_0);
    _rtB->B_97_1250_0 = _rtP->P_1395 * rtb_B_97_1173_0;
}
i = ssIsSampleHit(S, 1, 0);
if ((i != 0) && (_rtDW->RateTransition1_semaphoreTaken == 0)) {
    _rtDW->RateTransition1_Buffer0 = _rtB->B_97_1244_0;
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtDW->RateTransition1_semaphoreTaken = 1;
    rtb_B_97_1_0 = _rtDW->RateTransition1_Buffer0;
    _rtDW->RateTransition1_semaphoreTaken = 0;
    rtb_B_97_1_0 *= _rtP->P_1397;
    _rtB->B_97_1257_0 = ((rtb_B_97_1_0 / (rtb_B_97_17_0 + 2.2204460492503131e-16) - rtb_B_97_1236_0) - _rtP->P_1398 * rtb_B_97_17_0) * _rtP->P_1399;
    _rtB->B_97_1258_0 = _rtP->P_1400 * rtb_B_97_16_0;
    _rtB->B_97_1259_0 = _rtP->P_1401 * rtb_B_97_17_0;
    _rtB->B_97_295_0[0] = _rtP->P_1402 * _rtB->B_97_116_0[39] * _rtP->P_1405;
    _rtB->B_97_295_0[1] = _rtP->P_1403 * _rtB->B_97_116_0[40] * _rtP->P_1405;
    _rtB->B_97_295_0[2] = _rtP->P_1404 * _rtB->B_97_116_0[41] * _rtP->P_1405;
    if (_rtDW->systemEnable_n != 0) {
        _rtDW->lastSin_b = muDoubleScalarSin(_rtP->P_1408 * ssGetTaskTime(S, 2));
        _rtDW->lastCos_p = muDoubleScalarCos(_rtP->P_1408 * ssGetTaskTime(S, 2));
        _rtDW->systemEnable_n = 0;
    }
    rtb_B_97_16_0 = ((_rtDW->lastSin_b * _rtP->P_1412 + _rtDW->lastCos_p * _rtP->P_1411) * _rtP->P_1410 + (_rtDW->lastCos_p * _rtP->P_1412 - _rtDW->lastSin_b * _rtP->P_1411) * _rtP->P_1409) * _rtP->P_1406 + _rtP->P_1407;
    _rtB->B_97_1270_0[0] = _rtB->B_97_295_0[0] * rtb_B_97_16_0;
    _rtB->B_97_1270_0[1] = _rtB->B_97_295_0[1] * rtb_B_97_16_0;
    _rtB->B_97_1270_0[2] = _rtB->B_97_295_0[2] * rtb_B_97_16_0;
    if (_rtDW->Integ4_SYSTEM_ENABLE_c2 != 0) {
        _rtB->B_97_1271_0[0] = _rtDW->Integ4_DSTATE_dr[0];
        _rtB->B_97_1271_0[1] = _rtDW->Integ4_DSTATE_dr[1];
        _rtB->B_97_1271_0[2] = _rtDW->Integ4_DSTATE_dr[2];
    } else {
        _rtB->B_97_1271_0[0] = _rtP->P_1413 * _rtB->B_97_1270_0[0] + _rtDW->Integ4_DSTATE_dr[0];
        _rtB->B_97_1271_0[1] = _rtP->P_1413 * _rtB->B_97_1270_0[1] + _rtDW->Integ4_DSTATE_dr[1];
        _rtB->B_97_1271_0[2] = _rtP->P_1413 * _rtB->B_97_1270_0[2] + _rtDW->Integ4_DSTATE_dr[2];
    }
    _rtB->B_97_1272_0 = _rtP->P_1415;
                
          /* Level2 S-Function Block: '<S736>/B_82_239' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 1273, SS_CALL_MDL_OUTPUTS);

   
  
    if (ssGetTaskTime(S, 2) >= _rtB->B_97_1272_0) {
        _rtB->B_97_1280_0[0] = (_rtB->B_97_1271_0[0] - _rtB->B_97_1273_0[0]) * _rtP->P_1424 + (_rtP->P_412 * _rtB->B_97_1270_0[0] - _rtP->P_411 * _rtDW->UnitDelay_DSTATE_is[0]);
        _rtB->B_97_1280_0[1] = (_rtB->B_97_1271_0[1] - _rtB->B_97_1273_0[1]) * _rtP->P_1424 + (_rtP->P_412 * _rtB->B_97_1270_0[1] - _rtP->P_411 * _rtDW->UnitDelay_DSTATE_is[1]);
        _rtB->B_97_1280_0[2] = (_rtB->B_97_1271_0[2] - _rtB->B_97_1273_0[2]) * _rtP->P_1424 + (_rtP->P_412 * _rtB->B_97_1270_0[2] - _rtP->P_411 * _rtDW->UnitDelay_DSTATE_is[2]);
    } else {
        _rtB->B_97_1280_0[0] = _rtDW->UnitDelay1_DSTATE_cp[0];
        _rtB->B_97_1280_0[1] = _rtDW->UnitDelay1_DSTATE_cp[1];
        _rtB->B_97_1280_0[2] = _rtDW->UnitDelay1_DSTATE_cp[2];
    }
    if (_rtDW->systemEnable_h != 0) {
        _rtDW->lastSin_jg = muDoubleScalarSin(_rtP->P_1429 * ssGetTaskTime(S, 2));
        _rtDW->lastCos_fb = muDoubleScalarCos(_rtP->P_1429 * ssGetTaskTime(S, 2));
        _rtDW->systemEnable_h = 0;
    }
    rtb_B_97_16_0 = ((_rtDW->lastSin_jg * _rtP->P_1433 + _rtDW->lastCos_fb * _rtP->P_1432) * _rtP->P_1431 + (_rtDW->lastCos_fb * _rtP->P_1433 - _rtDW->lastSin_jg * _rtP->P_1432) * _rtP->P_1430) * _rtP->P_1427 + _rtP->P_1428;
    _rtB->B_97_1282_0[0] = _rtB->B_97_295_0[0] * rtb_B_97_16_0;
    _rtB->B_97_1282_0[1] = _rtB->B_97_295_0[1] * rtb_B_97_16_0;
    _rtB->B_97_1282_0[2] = _rtB->B_97_295_0[2] * rtb_B_97_16_0;
    if (_rtDW->Integ4_SYSTEM_ENABLE_jh != 0) {
        _rtB->B_97_1283_0[0] = _rtDW->Integ4_DSTATE_g5[0];
        _rtB->B_97_1283_0[1] = _rtDW->Integ4_DSTATE_g5[1];
        _rtB->B_97_1283_0[2] = _rtDW->Integ4_DSTATE_g5[2];
    } else {
        _rtB->B_97_1283_0[0] = _rtP->P_1434 * _rtB->B_97_1282_0[0] + _rtDW->Integ4_DSTATE_g5[0];
        _rtB->B_97_1283_0[1] = _rtP->P_1434 * _rtB->B_97_1282_0[1] + _rtDW->Integ4_DSTATE_g5[1];
        _rtB->B_97_1283_0[2] = _rtP->P_1434 * _rtB->B_97_1282_0[2] + _rtDW->Integ4_DSTATE_g5[2];
    }
    _rtB->B_97_1284_0 = _rtP->P_1436;
                
          /* Level2 S-Function Block: '<S734>/B_82_241' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 1285, SS_CALL_MDL_OUTPUTS);

   
  
    if (ssGetTaskTime(S, 2) >= _rtB->B_97_1284_0) {
        _rtB->B_97_1292_0[0] = (_rtB->B_97_1283_0[0] - _rtB->B_97_1285_0[0]) * _rtP->P_1445 + (_rtP->P_410 * _rtB->B_97_1282_0[0] - _rtP->P_409 * _rtDW->UnitDelay_DSTATE_ox[0]);
        _rtB->B_97_1292_0[1] = (_rtB->B_97_1283_0[1] - _rtB->B_97_1285_0[1]) * _rtP->P_1445 + (_rtP->P_410 * _rtB->B_97_1282_0[1] - _rtP->P_409 * _rtDW->UnitDelay_DSTATE_ox[1]);
        _rtB->B_97_1292_0[2] = (_rtB->B_97_1283_0[2] - _rtB->B_97_1285_0[2]) * _rtP->P_1445 + (_rtP->P_410 * _rtB->B_97_1282_0[2] - _rtP->P_409 * _rtDW->UnitDelay_DSTATE_ox[2]);
    } else {
        _rtB->B_97_1292_0[0] = _rtDW->UnitDelay1_DSTATE_f[0];
        _rtB->B_97_1292_0[1] = _rtDW->UnitDelay1_DSTATE_f[1];
        _rtB->B_97_1292_0[2] = _rtDW->UnitDelay1_DSTATE_f[2];
    }
    B_97_726_1_idx_0 = muDoubleScalarHypot(_rtB->B_97_1280_0[0], _rtB->B_97_1292_0[0]);
    B_97_719_0_idx_0 = muDoubleScalarAtan2(_rtB->B_97_1292_0[0], _rtB->B_97_1280_0[0]);
    B_97_726_1_idx_1 = muDoubleScalarHypot(_rtB->B_97_1280_0[1], _rtB->B_97_1292_0[1]);
    B_97_719_0_idx_1 = muDoubleScalarAtan2(_rtB->B_97_1292_0[1], _rtB->B_97_1280_0[1]);
    B_97_726_1_idx_2 = muDoubleScalarHypot(_rtB->B_97_1280_0[2], _rtB->B_97_1292_0[2]);
    B_97_719_0_idx_2 = muDoubleScalarAtan2(_rtB->B_97_1292_0[2], _rtB->B_97_1280_0[2]);
    _rtB->B_97_265_0[0] = _rtP->P_1448 * _rtB->B_97_116_0[59] * _rtP->P_1451;
    _rtB->B_97_265_0[1] = _rtP->P_1449 * _rtB->B_97_116_0[60] * _rtP->P_1451;
    _rtB->B_97_265_0[2] = _rtP->P_1450 * _rtB->B_97_116_0[61] * _rtP->P_1451;
    if (_rtDW->systemEnable_az != 0) {
        _rtDW->lastSin_g = muDoubleScalarSin(_rtP->P_1454 * ssGetTaskTime(S, 2));
        _rtDW->lastCos_b = muDoubleScalarCos(_rtP->P_1454 * ssGetTaskTime(S, 2));
        _rtDW->systemEnable_az = 0;
    }
    rtb_B_97_16_0 = ((_rtDW->lastSin_g * _rtP->P_1458 + _rtDW->lastCos_b * _rtP->P_1457) * _rtP->P_1456 + (_rtDW->lastCos_b * _rtP->P_1458 - _rtDW->lastSin_g * _rtP->P_1457) * _rtP->P_1455) * _rtP->P_1452 + _rtP->P_1453;
    _rtB->B_97_1300_0[0] = _rtB->B_97_265_0[0] * rtb_B_97_16_0;
    _rtB->B_97_1300_0[1] = _rtB->B_97_265_0[1] * rtb_B_97_16_0;
    _rtB->B_97_1300_0[2] = _rtB->B_97_265_0[2] * rtb_B_97_16_0;
    if (_rtDW->Integ4_SYSTEM_ENABLE_n != 0) {
        _rtB->B_97_1301_0[0] = _rtDW->Integ4_DSTATE_h3[0];
        _rtB->B_97_1301_0[1] = _rtDW->Integ4_DSTATE_h3[1];
        _rtB->B_97_1301_0[2] = _rtDW->Integ4_DSTATE_h3[2];
    } else {
        _rtB->B_97_1301_0[0] = _rtP->P_1459 * _rtB->B_97_1300_0[0] + _rtDW->Integ4_DSTATE_h3[0];
        _rtB->B_97_1301_0[1] = _rtP->P_1459 * _rtB->B_97_1300_0[1] + _rtDW->Integ4_DSTATE_h3[1];
        _rtB->B_97_1301_0[2] = _rtP->P_1459 * _rtB->B_97_1300_0[2] + _rtDW->Integ4_DSTATE_h3[2];
    }
    _rtB->B_97_1302_0 = _rtP->P_1461;
                
          /* Level2 S-Function Block: '<S742>/B_82_243' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 1303, SS_CALL_MDL_OUTPUTS);

   
  
    if (ssGetTaskTime(S, 2) >= _rtB->B_97_1302_0) {
        _rtB->B_97_1310_0[0] = (_rtB->B_97_1301_0[0] - _rtB->B_97_1303_0[0]) * _rtP->P_1470 + (_rtP->P_416 * _rtB->B_97_1300_0[0] - _rtP->P_415 * _rtDW->UnitDelay_DSTATE_bp[0]);
        _rtB->B_97_1310_0[1] = (_rtB->B_97_1301_0[1] - _rtB->B_97_1303_0[1]) * _rtP->P_1470 + (_rtP->P_416 * _rtB->B_97_1300_0[1] - _rtP->P_415 * _rtDW->UnitDelay_DSTATE_bp[1]);
        _rtB->B_97_1310_0[2] = (_rtB->B_97_1301_0[2] - _rtB->B_97_1303_0[2]) * _rtP->P_1470 + (_rtP->P_416 * _rtB->B_97_1300_0[2] - _rtP->P_415 * _rtDW->UnitDelay_DSTATE_bp[2]);
    } else {
        _rtB->B_97_1310_0[0] = _rtDW->UnitDelay1_DSTATE_gq[0];
        _rtB->B_97_1310_0[1] = _rtDW->UnitDelay1_DSTATE_gq[1];
        _rtB->B_97_1310_0[2] = _rtDW->UnitDelay1_DSTATE_gq[2];
    }
    if (_rtDW->systemEnable_c != 0) {
        _rtDW->lastSin_i = muDoubleScalarSin(_rtP->P_1475 * ssGetTaskTime(S, 2));
        _rtDW->lastCos_e = muDoubleScalarCos(_rtP->P_1475 * ssGetTaskTime(S, 2));
        _rtDW->systemEnable_c = 0;
    }
    rtb_B_97_16_0 = ((_rtDW->lastSin_i * _rtP->P_1479 + _rtDW->lastCos_e * _rtP->P_1478) * _rtP->P_1477 + (_rtDW->lastCos_e * _rtP->P_1479 - _rtDW->lastSin_i * _rtP->P_1478) * _rtP->P_1476) * _rtP->P_1473 + _rtP->P_1474;
    _rtB->B_97_1312_0[0] = _rtB->B_97_265_0[0] * rtb_B_97_16_0;
    _rtB->B_97_1312_0[1] = _rtB->B_97_265_0[1] * rtb_B_97_16_0;
    _rtB->B_97_1312_0[2] = _rtB->B_97_265_0[2] * rtb_B_97_16_0;
    if (_rtDW->Integ4_SYSTEM_ENABLE_cn != 0) {
        _rtB->B_97_1313_0[0] = _rtDW->Integ4_DSTATE_bj[0];
        _rtB->B_97_1313_0[1] = _rtDW->Integ4_DSTATE_bj[1];
        _rtB->B_97_1313_0[2] = _rtDW->Integ4_DSTATE_bj[2];
    } else {
        _rtB->B_97_1313_0[0] = _rtP->P_1480 * _rtB->B_97_1312_0[0] + _rtDW->Integ4_DSTATE_bj[0];
        _rtB->B_97_1313_0[1] = _rtP->P_1480 * _rtB->B_97_1312_0[1] + _rtDW->Integ4_DSTATE_bj[1];
        _rtB->B_97_1313_0[2] = _rtP->P_1480 * _rtB->B_97_1312_0[2] + _rtDW->Integ4_DSTATE_bj[2];
    }
    _rtB->B_97_1314_0 = _rtP->P_1482;
                
          /* Level2 S-Function Block: '<S740>/B_82_245' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 1315, SS_CALL_MDL_OUTPUTS);

   
  
    if (ssGetTaskTime(S, 2) >= _rtB->B_97_1314_0) {
        _rtB->B_97_1322_0[0] = (_rtB->B_97_1313_0[0] - _rtB->B_97_1315_0[0]) * _rtP->P_1491 + (_rtP->P_414 * _rtB->B_97_1312_0[0] - _rtP->P_413 * _rtDW->UnitDelay_DSTATE_po[0]);
        _rtB->B_97_1322_0[1] = (_rtB->B_97_1313_0[1] - _rtB->B_97_1315_0[1]) * _rtP->P_1491 + (_rtP->P_414 * _rtB->B_97_1312_0[1] - _rtP->P_413 * _rtDW->UnitDelay_DSTATE_po[1]);
        _rtB->B_97_1322_0[2] = (_rtB->B_97_1313_0[2] - _rtB->B_97_1315_0[2]) * _rtP->P_1491 + (_rtP->P_414 * _rtB->B_97_1312_0[2] - _rtP->P_413 * _rtDW->UnitDelay_DSTATE_po[2]);
    } else {
        _rtB->B_97_1322_0[0] = _rtDW->UnitDelay1_DSTATE_gu[0];
        _rtB->B_97_1322_0[1] = _rtDW->UnitDelay1_DSTATE_gu[1];
        _rtB->B_97_1322_0[2] = _rtDW->UnitDelay1_DSTATE_gu[2];
    }
    B_97_726_1_idx_0 = B_97_726_1_idx_0 * muDoubleScalarHypot(_rtB->B_97_1310_0[0], _rtB->B_97_1322_0[0]) * _rtP->P_1494;
    B_97_719_0_idx_0 = (_rtP->P_1495 * B_97_719_0_idx_0 - _rtP->P_1496 * muDoubleScalarAtan2(_rtB->B_97_1322_0[0], _rtB->B_97_1310_0[0])) * _rtP->P_1497;
    B_97_726_1_idx_1 = B_97_726_1_idx_1 * muDoubleScalarHypot(_rtB->B_97_1310_0[1], _rtB->B_97_1322_0[1]) * _rtP->P_1494;
    B_97_719_0_idx_1 = (_rtP->P_1495 * B_97_719_0_idx_1 - _rtP->P_1496 * muDoubleScalarAtan2(_rtB->B_97_1322_0[1], _rtB->B_97_1310_0[1])) * _rtP->P_1497;
    B_97_726_1_idx_2 = B_97_726_1_idx_2 * muDoubleScalarHypot(_rtB->B_97_1310_0[2], _rtB->B_97_1322_0[2]) * _rtP->P_1494;
    rtb_B_97_1_0 = (_rtP->P_1495 * B_97_719_0_idx_2 - _rtP->P_1496 * muDoubleScalarAtan2(_rtB->B_97_1322_0[2], _rtB->B_97_1310_0[2])) * _rtP->P_1497;
    muDoubleScalarSinCos(B_97_719_0_idx_0, &B_97_1331_0_idx_0, &rtb_B_97_1173_0);
    muDoubleScalarSinCos(B_97_719_0_idx_1, &B_97_1331_0_idx_1, &rtb_B_97_1236_0);
    muDoubleScalarSinCos(rtb_B_97_1_0, &B_97_1331_0_idx_2, &rtb_B_97_5_0);
    _rtB->B_97_1333_0 = (B_97_726_1_idx_0 * rtb_B_97_1173_0 + B_97_726_1_idx_1 * rtb_B_97_1236_0) + B_97_726_1_idx_2 * rtb_B_97_5_0;
}
_rtB->B_97_1334_0 = _rtX->integ1_CSTATE_oz;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T_PWORK_g.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T_PWORK_g.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_1499;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        B_97_1335_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T_IWORK_b.CircularBufSize,
        &_rtDW->T_IWORK_b.Last,
        _rtDW->T_IWORK_b.Tail,
        _rtDW->T_IWORK_b.Head,
        _rtP->P_1500,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
B_97_1340_0.re = _rtB->B_97_1334_0 - B_97_1335_0;
_rtB->B_97_1337_0 = _rtX->Integ2_CSTATE_k;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK_bd.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK_bd.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_1502;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        B_97_1338_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T1_IWORK_fg.CircularBufSize,
        &_rtDW->T1_IWORK_fg.Last,
        _rtDW->T1_IWORK_fg.Tail,
        _rtDW->T1_IWORK_fg.Head,
        _rtP->P_1503,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
B_97_1340_0.im = _rtB->B_97_1337_0 - B_97_1338_0;
_rtB->B_97_1341_0 = _rtX->integ1_CSTATE_nb;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T_PWORK_d.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T_PWORK_d.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_1505;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        B_97_1342_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T_IWORK_m.CircularBufSize,
        &_rtDW->T_IWORK_m.Last,
        _rtDW->T_IWORK_m.Tail,
        _rtDW->T_IWORK_m.Head,
        _rtP->P_1506,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
B_97_1347_0.re = _rtB->B_97_1341_0 - B_97_1342_0;
_rtB->B_97_1344_0 = _rtX->Integ2_CSTATE_ik;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK_bk.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK_bk.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_1508;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        B_97_1345_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T1_IWORK_m.CircularBufSize,
        &_rtDW->T1_IWORK_m.Last,
        _rtDW->T1_IWORK_m.Tail,
        _rtDW->T1_IWORK_m.Head,
        _rtP->P_1509,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
B_97_1347_0.im = _rtB->B_97_1344_0 - B_97_1345_0;
_rtB->B_97_1348_0 = _rtX->integ1_CSTATE_g;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T_PWORK_l4.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T_PWORK_l4.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_1511;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        B_97_1349_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T_IWORK_c1.CircularBufSize,
        &_rtDW->T_IWORK_c1.Last,
        _rtDW->T_IWORK_c1.Tail,
        _rtDW->T_IWORK_c1.Head,
        _rtP->P_1512,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
B_97_1354_0.re = _rtB->B_97_1348_0 - B_97_1349_0;
_rtB->B_97_1351_0 = _rtX->Integ2_CSTATE_o;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK_p.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK_p.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_1514;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        B_97_1352_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T1_IWORK_c4.CircularBufSize,
        &_rtDW->T1_IWORK_c4.Last,
        _rtDW->T1_IWORK_c4.Tail,
        _rtDW->T1_IWORK_c4.Head,
        _rtP->P_1515,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
B_97_1354_0.im = _rtB->B_97_1351_0 - B_97_1352_0;
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    _rtB->B_97_1356_0 = _rtB->B_97_1355_0;
}





      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_NegSeqComputation(S, _rtB->B_97_1356_0, B_97_1340_0, B_97_1347_0, B_97_1354_0, &_rtB->PosSeqComputation_b, &_rtDW->PosSeqComputation_b, &_rtP->PosSeqComputation_b);
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    _rtB->B_97_1359_0 = _rtB->B_97_1358_0;
}





      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_NegSeqComputation(S, _rtB->B_97_1359_0, B_97_1340_0, B_97_1347_0, B_97_1354_0, &_rtB->NegSeqComputation_g, &_rtDW->NegSeqComputation_g, &_rtP->NegSeqComputation_g);
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    _rtB->B_97_1362_0 = _rtB->B_97_1361_0;
}





      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_ZeroSeqComputation(S, _rtB->B_97_1362_0, B_97_1340_0, B_97_1347_0, B_97_1354_0, &_rtB->ZeroSeqComputation_n, &_rtDW->ZeroSeqComputation_n, &_rtP->ZeroSeqComputation_n);
_rtB->B_97_1364_0[0] = muDoubleScalarHypot(_rtB->PosSeqComputation_b.B_39_2_0.re, _rtB->PosSeqComputation_b.B_39_2_0.im);
_rtB->B_97_1364_0[1] = muDoubleScalarHypot(_rtB->NegSeqComputation_g.B_39_2_0.re, _rtB->NegSeqComputation_g.B_39_2_0.im);
_rtB->B_97_1364_0[2] = muDoubleScalarHypot(_rtB->ZeroSeqComputation_n.B_41_1_0.re, _rtB->ZeroSeqComputation_n.B_41_1_0.im);
rtb_B_97_30_0[0] = muDoubleScalarAtan2(_rtB->PosSeqComputation_b.B_39_2_0.im, _rtB->PosSeqComputation_b.B_39_2_0.re);
rtb_B_97_30_0[1] = muDoubleScalarAtan2(_rtB->NegSeqComputation_g.B_39_2_0.im, _rtB->NegSeqComputation_g.B_39_2_0.re);
rtb_B_97_30_0[2] = muDoubleScalarAtan2(_rtB->ZeroSeqComputation_n.B_41_1_0.im, _rtB->ZeroSeqComputation_n.B_41_1_0.re);
_rtB->B_97_1365_0 = _rtP->P_1519 * _rtB->B_97_1364_0[0];
_rtB->B_97_1366_0 = 0.0;
_rtB->B_97_1366_0 += _rtP->P_1521 * _rtX->TransferFcn1_CSTATE_f;
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    B_97_726_1_idx_0 *= B_97_1331_0_idx_0;
    B_97_726_1_idx_1 *= B_97_1331_0_idx_1;
    B_97_726_1_idx_2 *= B_97_1331_0_idx_2;
    _rtB->B_97_1368_0 = (B_97_726_1_idx_0 + B_97_726_1_idx_1) + B_97_726_1_idx_2;
}
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 1369, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 1370, SS_CALL_MDL_OUTPUTS);
   
  
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 1371, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 1372, SS_CALL_MDL_OUTPUTS);
   
  
}
_rtB->B_97_1374_0 = _rtP->P_1523 * _rtB->B_97_1364_0[0];
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    if (ssIsMajorTimeStep(S) != 0) {
        _rtDW->RelationalOperator_Mode_ej[0] = (_rtB->B_97_1364_0[1] > _rtB->B_97_1374_0);
        _rtDW->RelationalOperator_Mode_ej[1] = (_rtB->B_97_1364_0[2] > _rtB->B_97_1374_0);
    }
    _rtB->B_97_1375_0[0] = _rtDW->RelationalOperator_Mode_ej[0];
    _rtB->B_97_1375_0[1] = _rtDW->RelationalOperator_Mode_ej[1];
}
if (_rtB->B_97_1375_0[0] >= _rtP->P_1524) {
    B_97_1376_0[0] = rtb_B_97_30_0[1];
} else {
    B_97_1376_0[0] = _rtB->B_97_1373_0;
}
if (_rtB->B_97_1375_0[1] >= _rtP->P_1524) {
    B_97_1376_0[1] = rtb_B_97_30_0[2];
} else {
    B_97_1376_0[1] = _rtB->B_97_1373_0;
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_1377_0 = _rtP->P_1525 * _rtB->B_97_295_0[0];
}
rtb_B_97_16_0 = muDoubleScalarSin(_rtP->P_1528 * ssGetTaskTime(S, 0) + _rtP->P_1529) * _rtP->P_1526 + _rtP->P_1527;
_rtB->B_97_1379_0 = rtb_B_97_16_0 * _rtB->B_97_1377_0;
rtb_B_97_17_0 = muDoubleScalarSin(_rtP->P_1532 * ssGetTaskTime(S, 0) + _rtP->P_1533) * _rtP->P_1530 + _rtP->P_1531;
_rtB->B_97_1381_0 = _rtB->B_97_1377_0 * rtb_B_97_17_0;
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_1382_0 = _rtP->P_1534 * _rtB->B_97_295_0[1];
}
_rtB->B_97_1383_0 = rtb_B_97_16_0 * _rtB->B_97_1382_0;
_rtB->B_97_1384_0 = _rtB->B_97_1382_0 * rtb_B_97_17_0;
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_1385_0 = _rtP->P_1535 * _rtB->B_97_295_0[2];
}
_rtB->B_97_1386_0 = rtb_B_97_16_0 * _rtB->B_97_1385_0;
_rtB->B_97_1387_0 = _rtB->B_97_1385_0 * rtb_B_97_17_0;
B_97_1388_0 = _rtP->P_1536 * rtb_B_97_30_0[0];
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    B_97_321_0_idx_0 = _rtP->P_1537 * _rtB->B_97_295_0[0];
    B_97_321_0_idx_1 = _rtP->P_1537 * _rtB->B_97_295_0[1];
    B_97_321_0_idx_2 = _rtP->P_1537 * _rtB->B_97_295_0[2];
    rtb_B_97_17_0 = muDoubleScalarMod(_rtDW->DiscreteTimeIntegrator_DSTATE_fi, _rtB->B_97_1391_0);
    rtb_B_97_16_0 = muDoubleScalarSin(rtb_B_97_17_0);
    rtb_B_97_17_0 = muDoubleScalarCos(rtb_B_97_17_0);
    rtb_B_97_1236_0 = (0.0 - rtb_B_97_16_0 * _rtB->B_97_1395_0) - rtb_B_97_17_0 * _rtB->B_97_1397_0;
    rtb_B_97_1_0 = rtb_B_97_16_0 * _rtB->B_97_1397_0 - rtb_B_97_17_0 * _rtB->B_97_1395_0;
    rtb_B_97_1173_0 = (0.0 - rtb_B_97_1236_0) - rtb_B_97_16_0;
    rtb_B_97_5_0 = (0.0 - rtb_B_97_1_0) - rtb_B_97_17_0;
    _rtB->B_97_1405_0 = ((B_97_321_0_idx_0 * rtb_B_97_17_0 + B_97_321_0_idx_1 * rtb_B_97_1_0) + B_97_321_0_idx_2 * rtb_B_97_5_0) * 0.66666666666666663;
    if (_rtDW->Integ4_SYSTEM_ENABLE_kv != 0) {
        _rtB->B_97_1406_0 = _rtDW->Integ4_DSTATE_f2;
    } else {
        _rtB->B_97_1406_0 = _rtP->P_1543 * _rtB->B_97_1405_0 + _rtDW->Integ4_DSTATE_f2;
    }
    _rtB->B_97_1407_0 = _rtDW->UnitDelay_DSTATE_l3;
    if (_rtB->B_97_1407_0 > _rtP->P_1546) {
        B_5_0_0_idx_3 = _rtP->P_1546;
    } else if (_rtB->B_97_1407_0 < _rtP->P_1547) {
        B_5_0_0_idx_3 = _rtP->P_1547;
    } else {
        B_5_0_0_idx_3 = _rtB->B_97_1407_0;
    }
    rtb_B_97_5_1 = 1.0 / B_5_0_0_idx_3 / 1.0e-5;
    rtb_B_97_27_0 = muDoubleScalarCeil(rtb_B_97_5_1);
    _rtB->B_97_1411_0 = _rtP->P_1548 * rtb_B_97_27_0;
                
          /* Level2 S-Function Block: '<S728>/B_82_260' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 1412, SS_CALL_MDL_OUTPUTS);

   
  
    rtb_B_97_82_0 = ssGetTaskTime(S, 2);
    if (rtb_B_97_82_0 < _rtP->P_1550) {
        rtb_B_97_82_0 = _rtP->P_1551;
    } else {
        rtb_B_97_82_0 = _rtP->P_1552;
    }
    if (rtb_B_97_82_0 >= _rtP->P_1554) {
        rtb_B_97_27_0 = rtb_B_97_5_1 - rtb_B_97_27_0;
        rtb_B_97_5_1 = ((_rtB->B_97_1405_0 - _rtDW->UnitDelay_DSTATE_gl) * rtb_B_97_27_0 * _rtP->P_408 + _rtB->B_97_1405_0) * (rtb_B_97_27_0 / rtb_B_97_5_1) + (_rtB->B_97_1406_0 - _rtB->B_97_1412_0) * _rtB->B_97_1407_0;
    } else {
        rtb_B_97_5_1 = _rtP->P_1553;
    }
    _rtB->B_97_1419_0 = _rtB->B_97_1418_0;
    i = ssIsSampleHit(S, 2, 0);
    if ((i != 0) && (ssIsMajorTimeStep(S) != 0)) {
        if (_rtB->B_97_1419_0 > 0.0) {
            if (!_rtDW->AutomaticGainControl_MODE_a) {
                if (ssGetTaskTime(S, 2) != ssGetTStart(S)) {
                    ssSetBlockStateForSolverChangedAtMajorStep(S);
                }
  
                _rtDW->Integ4_SYSTEM_ENABLE_fs = 1U;
                _rtDW->Integ4_SYSTEM_ENABLE_mo = 1U;
                _rtDW->AutomaticGainControl_MODE_a = true;
            }
        } else {
            if (_rtDW->AutomaticGainControl_MODE_a) {
                ssSetBlockStateForSolverChangedAtMajorStep(S);
  
                _rtDW->AutomaticGainControl_MODE_a = false;
            }
        }
    }
    if (_rtDW->AutomaticGainControl_MODE_a) {
        i = ssIsSampleHit(S, 2, 0);
        if (i != 0) {
            rtb_B_97_27_0 = (0.0 - rtb_B_97_16_0 * _rtB->B_55_0_0) - rtb_B_97_17_0 * _rtB->B_55_2_0;
            rtb_B_97_28_0 = rtb_B_97_16_0 * _rtB->B_55_2_0 - rtb_B_97_17_0 * _rtB->B_55_0_0;
            B_97_726_1_idx_0 = (0.0 - rtb_B_97_27_0) - rtb_B_97_16_0;
            B_97_726_1_idx_1 = (0.0 - rtb_B_97_28_0) - rtb_B_97_17_0;
            _rtB->B_55_10_0 = ((B_97_321_0_idx_0 * rtb_B_97_16_0 + B_97_321_0_idx_1 * rtb_B_97_27_0) + B_97_321_0_idx_2 * B_97_726_1_idx_0) * 0.66666666666666663;
            if (_rtDW->Integ4_SYSTEM_ENABLE_fs != 0) {
                _rtB->B_55_11_0 = _rtDW->Integ4_DSTATE_i0;
            } else {
                _rtB->B_55_11_0 = _rtP->P_382 * _rtB->B_55_10_0 + _rtDW->Integ4_DSTATE_i0;
            }
            if (_rtB->B_97_1407_0 > _rtP->P_384) {
                B_5_0_0_idx_3 = _rtP->P_384;
            } else if (_rtB->B_97_1407_0 < _rtP->P_385) {
                B_5_0_0_idx_3 = _rtP->P_385;
            } else {
                B_5_0_0_idx_3 = _rtB->B_97_1407_0;
            }
            rtb_B_97_966_0 = 1.0 / B_5_0_0_idx_3 / 1.0e-5;
            rtb_B_97_967_0 = muDoubleScalarCeil(rtb_B_97_966_0);
            _rtB->B_55_15_0 = _rtP->P_386 * rtb_B_97_967_0;
                
          /* Level2 S-Function Block: '<S723>/B_50_0' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 55, 16, SS_CALL_MDL_OUTPUTS);

   
  
            rtb_B_97_82_0 = ssGetTaskTime(S, 2);
            if (rtb_B_97_82_0 < _rtP->P_388) {
                rtb_B_97_82_0 = _rtP->P_389;
            } else {
                rtb_B_97_82_0 = _rtP->P_390;
            }
            if (rtb_B_97_82_0 >= _rtP->P_392) {
                rtb_B_97_967_0 = rtb_B_97_966_0 - rtb_B_97_967_0;
                rtb_B_97_966_0 = ((_rtB->B_55_10_0 - _rtDW->UnitDelay_DSTATE_pv) * rtb_B_97_967_0 * _rtP->P_377 + _rtB->B_55_10_0) * (rtb_B_97_967_0 / rtb_B_97_966_0) + (_rtB->B_55_11_0 - _rtB->B_55_16_0) * _rtB->B_97_1407_0;
            } else {
                rtb_B_97_966_0 = _rtP->P_391;
            }
            _rtB->B_55_22_0 = ((B_97_321_0_idx_0 * rtb_B_97_17_0 + B_97_321_0_idx_1 * rtb_B_97_28_0) + B_97_321_0_idx_2 * B_97_726_1_idx_1) * 0.66666666666666663;
            if (_rtDW->Integ4_SYSTEM_ENABLE_mo != 0) {
                _rtB->B_55_23_0 = _rtDW->Integ4_DSTATE_c;
            } else {
                _rtB->B_55_23_0 = _rtP->P_393 * _rtB->B_55_22_0 + _rtDW->Integ4_DSTATE_c;
            }
            if (_rtB->B_97_1407_0 > _rtP->P_395) {
                B_5_0_0_idx_3 = _rtP->P_395;
            } else if (_rtB->B_97_1407_0 < _rtP->P_396) {
                B_5_0_0_idx_3 = _rtP->P_396;
            } else {
                B_5_0_0_idx_3 = _rtB->B_97_1407_0;
            }
            rtb_B_97_967_0 = 1.0 / B_5_0_0_idx_3 / 1.0e-5;
            rtb_B_97_972_0 = muDoubleScalarCeil(rtb_B_97_967_0);
            _rtB->B_55_27_0 = _rtP->P_397 * rtb_B_97_972_0;
                
          /* Level2 S-Function Block: '<S725>/B_50_2' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 55, 28, SS_CALL_MDL_OUTPUTS);

   
  
            rtb_B_97_82_0 = ssGetTaskTime(S, 2);
            if (rtb_B_97_82_0 < _rtP->P_399) {
                rtb_B_97_82_0 = _rtP->P_400;
            } else {
                rtb_B_97_82_0 = _rtP->P_401;
            }
            if (rtb_B_97_82_0 >= _rtP->P_403) {
                rtb_B_97_972_0 = rtb_B_97_967_0 - rtb_B_97_972_0;
                rtb_B_97_967_0 = ((_rtB->B_55_22_0 - _rtDW->UnitDelay_DSTATE_as) * rtb_B_97_972_0 * _rtP->P_378 + _rtB->B_55_22_0) * (rtb_B_97_972_0 / rtb_B_97_967_0) + (_rtB->B_55_23_0 - _rtB->B_55_28_0) * _rtB->B_97_1407_0;
            } else {
                rtb_B_97_967_0 = _rtP->P_402;
            }
            B_55_36_0 = _rtP->P_404 * muDoubleScalarAtan2(rtb_B_97_967_0, rtb_B_97_966_0);
            tmpForInput[0] = B_97_321_0_idx_0;
            tmpForInput[1] = B_97_321_0_idx_1;
            tmpForInput[2] = B_97_321_0_idx_2;
            tmpForInput[3] = rtb_B_97_16_0;
            tmpForInput[4] = rtb_B_97_17_0;
            tmpForInput[5] = rtb_B_97_27_0;
            tmpForInput[6] = rtb_B_97_28_0;
            tmpForInput[7] = B_97_726_1_idx_0;
            tmpForInput[8] = B_97_726_1_idx_1;
            B_97_1331_0_idx_0 = -0.0;
            for (i = 0; i < 9; i++) {
                B_97_1331_0_idx_0 += tmpForInput[i];
            }
            B_55_38_0 = _rtP->P_405 * B_97_1331_0_idx_0;
            B_5_0_0_idx_3 = muDoubleScalarHypot(rtb_B_97_966_0, rtb_B_97_967_0);
            if (B_5_0_0_idx_3 > _rtP->P_406) {
                B_5_0_0_idx_3 = _rtP->P_406;
            } else {
                if (B_5_0_0_idx_3 < _rtP->P_407) {
                    B_5_0_0_idx_3 = _rtP->P_407;
                }
            }
            _rtB->B_55_40_0 = 1.0 / B_5_0_0_idx_3;
        }
        if (ssIsMajorTimeStep(S) != 0) {
            srUpdateBC(_rtDW->AutomaticGainControl_SubsysRanBC_e);
        }
    }
    B_97_1331_0_idx_0 = rtb_B_97_5_1 * _rtB->B_55_40_0;
    _rtB->B_97_1426_0 = _rtP->P_1561 * B_97_1331_0_idx_0 * _rtP->P_1562;
    B_97_1427_0 = _rtDW->UD_DSTATE_aa;
    B_97_1428_0 = _rtB->B_97_1426_0 - B_97_1427_0;
    B_5_0_0_idx_3 = (_rtP->P_1556 * B_97_1331_0_idx_0 + _rtDW->DiscreteTimeIntegrator_DSTATE_a) + B_97_1428_0;
    if (B_5_0_0_idx_3 > _rtP->P_1564) {
        _rtB->B_97_1430_0 = _rtP->P_1564;
    } else if (B_5_0_0_idx_3 < _rtP->P_1565) {
        _rtB->B_97_1430_0 = _rtP->P_1565;
    } else {
        _rtB->B_97_1430_0 = B_5_0_0_idx_3;
    }
    B_5_0_0_idx_3 = _rtP->P_1566 * _rtB->B_97_1430_0 - _rtDW->UnitDelay_DSTATE_oi;
    if (B_5_0_0_idx_3 > _rtP->P_1568) {
        B_5_0_0_idx_3 = _rtP->P_1568;
    } else {
        if (B_5_0_0_idx_3 < _rtP->P_1569) {
            B_5_0_0_idx_3 = _rtP->P_1569;
        }
    }
    _rtB->B_97_1435_0 = B_5_0_0_idx_3 + _rtDW->UnitDelay_DSTATE_oi;
          {
                 _rtB->B_97_1436_0 = (_rtP->P_1572[0])*_rtDW->DiscreteStateSpace_DSTATE_a[0] + (_rtP->P_1572[1])*_rtDW->DiscreteStateSpace_DSTATE_a[1];
    _rtB->B_97_1436_0 += _rtP->P_1573*_rtB->B_97_1435_0;
  }

   
  
    _rtB->B_97_1438_0 = _rtP->P_1575 * B_97_1331_0_idx_0;
    B_97_1439_0 = ((B_97_321_0_idx_0 * rtb_B_97_16_0 + B_97_321_0_idx_1 * rtb_B_97_1236_0) + B_97_321_0_idx_2 * rtb_B_97_1173_0) * 0.66666666666666663;
    tmpForInput[0] = B_97_321_0_idx_0;
    tmpForInput[1] = B_97_321_0_idx_1;
    tmpForInput[2] = B_97_321_0_idx_2;
    tmpForInput[3] = rtb_B_97_16_0;
    tmpForInput[4] = rtb_B_97_17_0;
    tmpForInput[5] = rtb_B_97_1236_0;
    tmpForInput[6] = rtb_B_97_1_0;
    tmpForInput[7] = rtb_B_97_1173_0;
    tmpForInput[8] = rtb_B_97_5_0;
    rtb_B_97_967_0 = -0.0;
    for (i = 0; i < 9; i++) {
        rtb_B_97_967_0 += tmpForInput[i];
    }
    B_97_1441_0 = _rtP->P_1576 * rtb_B_97_967_0;
    _rtB->B_97_295_0[0] = _rtP->P_1577 * _rtB->B_97_116_0[42] * _rtP->P_1580;
    _rtB->B_97_295_0[1] = _rtP->P_1578 * _rtB->B_97_116_0[43] * _rtP->P_1580;
    _rtB->B_97_295_0[2] = _rtP->P_1579 * _rtB->B_97_116_0[44] * _rtP->P_1580;
    if (_rtDW->systemEnable_hq != 0) {
        _rtDW->lastSin_im = muDoubleScalarSin(_rtP->P_1583 * ssGetTaskTime(S, 2));
        _rtDW->lastCos_k = muDoubleScalarCos(_rtP->P_1583 * ssGetTaskTime(S, 2));
        _rtDW->systemEnable_hq = 0;
    }
    rtb_B_97_967_0 = ((_rtDW->lastSin_im * _rtP->P_1587 + _rtDW->lastCos_k * _rtP->P_1586) * _rtP->P_1585 + (_rtDW->lastCos_k * _rtP->P_1587 - _rtDW->lastSin_im * _rtP->P_1586) * _rtP->P_1584) * _rtP->P_1581 + _rtP->P_1582;
    _rtB->B_97_1481_0[0] = _rtB->B_97_295_0[0] * rtb_B_97_967_0;
    _rtB->B_97_1481_0[1] = _rtB->B_97_295_0[1] * rtb_B_97_967_0;
    _rtB->B_97_1481_0[2] = _rtB->B_97_295_0[2] * rtb_B_97_967_0;
    if (_rtDW->Integ4_SYSTEM_ENABLE_pe != 0) {
        _rtB->B_97_1482_0[0] = _rtDW->Integ4_DSTATE_gh[0];
        _rtB->B_97_1482_0[1] = _rtDW->Integ4_DSTATE_gh[1];
        _rtB->B_97_1482_0[2] = _rtDW->Integ4_DSTATE_gh[2];
    } else {
        _rtB->B_97_1482_0[0] = _rtP->P_1588 * _rtB->B_97_1481_0[0] + _rtDW->Integ4_DSTATE_gh[0];
        _rtB->B_97_1482_0[1] = _rtP->P_1588 * _rtB->B_97_1481_0[1] + _rtDW->Integ4_DSTATE_gh[1];
        _rtB->B_97_1482_0[2] = _rtP->P_1588 * _rtB->B_97_1481_0[2] + _rtDW->Integ4_DSTATE_gh[2];
    }
    _rtB->B_97_1483_0 = _rtP->P_1590;
                
          /* Level2 S-Function Block: '<S820>/B_82_298' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 1484, SS_CALL_MDL_OUTPUTS);

   
  
    if (ssGetTaskTime(S, 2) >= _rtB->B_97_1483_0) {
        _rtB->B_97_1491_0[0] = (_rtB->B_97_1482_0[0] - _rtB->B_97_1484_0[0]) * _rtP->P_1599 + (_rtP->P_492 * _rtB->B_97_1481_0[0] - _rtP->P_491 * _rtDW->UnitDelay_DSTATE_bw[0]);
        _rtB->B_97_1491_0[1] = (_rtB->B_97_1482_0[1] - _rtB->B_97_1484_0[1]) * _rtP->P_1599 + (_rtP->P_492 * _rtB->B_97_1481_0[1] - _rtP->P_491 * _rtDW->UnitDelay_DSTATE_bw[1]);
        _rtB->B_97_1491_0[2] = (_rtB->B_97_1482_0[2] - _rtB->B_97_1484_0[2]) * _rtP->P_1599 + (_rtP->P_492 * _rtB->B_97_1481_0[2] - _rtP->P_491 * _rtDW->UnitDelay_DSTATE_bw[2]);
    } else {
        _rtB->B_97_1491_0[0] = _rtDW->UnitDelay1_DSTATE_hp[0];
        _rtB->B_97_1491_0[1] = _rtDW->UnitDelay1_DSTATE_hp[1];
        _rtB->B_97_1491_0[2] = _rtDW->UnitDelay1_DSTATE_hp[2];
    }
    if (_rtDW->systemEnable_ci != 0) {
        _rtDW->lastSin_bt = muDoubleScalarSin(_rtP->P_1604 * ssGetTaskTime(S, 2));
        _rtDW->lastCos_m = muDoubleScalarCos(_rtP->P_1604 * ssGetTaskTime(S, 2));
        _rtDW->systemEnable_ci = 0;
    }
    rtb_B_97_967_0 = ((_rtDW->lastSin_bt * _rtP->P_1608 + _rtDW->lastCos_m * _rtP->P_1607) * _rtP->P_1606 + (_rtDW->lastCos_m * _rtP->P_1608 - _rtDW->lastSin_bt * _rtP->P_1607) * _rtP->P_1605) * _rtP->P_1602 + _rtP->P_1603;
    _rtB->B_97_1493_0[0] = _rtB->B_97_295_0[0] * rtb_B_97_967_0;
    _rtB->B_97_1493_0[1] = _rtB->B_97_295_0[1] * rtb_B_97_967_0;
    _rtB->B_97_1493_0[2] = _rtB->B_97_295_0[2] * rtb_B_97_967_0;
    if (_rtDW->Integ4_SYSTEM_ENABLE_ix != 0) {
        _rtB->B_97_1494_0[0] = _rtDW->Integ4_DSTATE_n4[0];
        _rtB->B_97_1494_0[1] = _rtDW->Integ4_DSTATE_n4[1];
        _rtB->B_97_1494_0[2] = _rtDW->Integ4_DSTATE_n4[2];
    } else {
        _rtB->B_97_1494_0[0] = _rtP->P_1609 * _rtB->B_97_1493_0[0] + _rtDW->Integ4_DSTATE_n4[0];
        _rtB->B_97_1494_0[1] = _rtP->P_1609 * _rtB->B_97_1493_0[1] + _rtDW->Integ4_DSTATE_n4[1];
        _rtB->B_97_1494_0[2] = _rtP->P_1609 * _rtB->B_97_1493_0[2] + _rtDW->Integ4_DSTATE_n4[2];
    }
    _rtB->B_97_1495_0 = _rtP->P_1611;
                
          /* Level2 S-Function Block: '<S818>/B_82_300' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 1496, SS_CALL_MDL_OUTPUTS);

   
  
    if (ssGetTaskTime(S, 2) >= _rtB->B_97_1495_0) {
        _rtB->B_97_1503_0[0] = (_rtB->B_97_1494_0[0] - _rtB->B_97_1496_0[0]) * _rtP->P_1620 + (_rtP->P_490 * _rtB->B_97_1493_0[0] - _rtP->P_489 * _rtDW->UnitDelay_DSTATE_kh[0]);
        _rtB->B_97_1503_0[1] = (_rtB->B_97_1494_0[1] - _rtB->B_97_1496_0[1]) * _rtP->P_1620 + (_rtP->P_490 * _rtB->B_97_1493_0[1] - _rtP->P_489 * _rtDW->UnitDelay_DSTATE_kh[1]);
        _rtB->B_97_1503_0[2] = (_rtB->B_97_1494_0[2] - _rtB->B_97_1496_0[2]) * _rtP->P_1620 + (_rtP->P_490 * _rtB->B_97_1493_0[2] - _rtP->P_489 * _rtDW->UnitDelay_DSTATE_kh[2]);
    } else {
        _rtB->B_97_1503_0[0] = _rtDW->UnitDelay1_DSTATE_fp[0];
        _rtB->B_97_1503_0[1] = _rtDW->UnitDelay1_DSTATE_fp[1];
        _rtB->B_97_1503_0[2] = _rtDW->UnitDelay1_DSTATE_fp[2];
    }
    B_97_726_1_idx_0 = muDoubleScalarHypot(_rtB->B_97_1491_0[0], _rtB->B_97_1503_0[0]);
    B_97_719_0_idx_0 = muDoubleScalarAtan2(_rtB->B_97_1503_0[0], _rtB->B_97_1491_0[0]);
    B_97_726_1_idx_1 = muDoubleScalarHypot(_rtB->B_97_1491_0[1], _rtB->B_97_1503_0[1]);
    B_97_719_0_idx_1 = muDoubleScalarAtan2(_rtB->B_97_1503_0[1], _rtB->B_97_1491_0[1]);
    B_97_726_1_idx_2 = muDoubleScalarHypot(_rtB->B_97_1491_0[2], _rtB->B_97_1503_0[2]);
    B_97_719_0_idx_2 = muDoubleScalarAtan2(_rtB->B_97_1503_0[2], _rtB->B_97_1491_0[2]);
    _rtB->B_97_265_0[0] = _rtP->P_1623 * _rtB->B_97_116_0[62] * _rtP->P_1626;
    _rtB->B_97_265_0[1] = _rtP->P_1624 * _rtB->B_97_116_0[63] * _rtP->P_1626;
    _rtB->B_97_265_0[2] = _rtP->P_1625 * _rtB->B_97_116_0[64] * _rtP->P_1626;
    if (_rtDW->systemEnable_hi != 0) {
        _rtDW->lastSin_p = muDoubleScalarSin(_rtP->P_1629 * ssGetTaskTime(S, 2));
        _rtDW->lastCos_fk = muDoubleScalarCos(_rtP->P_1629 * ssGetTaskTime(S, 2));
        _rtDW->systemEnable_hi = 0;
    }
    rtb_B_97_967_0 = ((_rtDW->lastSin_p * _rtP->P_1633 + _rtDW->lastCos_fk * _rtP->P_1632) * _rtP->P_1631 + (_rtDW->lastCos_fk * _rtP->P_1633 - _rtDW->lastSin_p * _rtP->P_1632) * _rtP->P_1630) * _rtP->P_1627 + _rtP->P_1628;
    _rtB->B_97_1511_0[0] = _rtB->B_97_265_0[0] * rtb_B_97_967_0;
    _rtB->B_97_1511_0[1] = _rtB->B_97_265_0[1] * rtb_B_97_967_0;
    _rtB->B_97_1511_0[2] = _rtB->B_97_265_0[2] * rtb_B_97_967_0;
    if (_rtDW->Integ4_SYSTEM_ENABLE_lx != 0) {
        _rtB->B_97_1512_0[0] = _rtDW->Integ4_DSTATE_pq[0];
        _rtB->B_97_1512_0[1] = _rtDW->Integ4_DSTATE_pq[1];
        _rtB->B_97_1512_0[2] = _rtDW->Integ4_DSTATE_pq[2];
    } else {
        _rtB->B_97_1512_0[0] = _rtP->P_1634 * _rtB->B_97_1511_0[0] + _rtDW->Integ4_DSTATE_pq[0];
        _rtB->B_97_1512_0[1] = _rtP->P_1634 * _rtB->B_97_1511_0[1] + _rtDW->Integ4_DSTATE_pq[1];
        _rtB->B_97_1512_0[2] = _rtP->P_1634 * _rtB->B_97_1511_0[2] + _rtDW->Integ4_DSTATE_pq[2];
    }
    _rtB->B_97_1513_0 = _rtP->P_1636;
                
          /* Level2 S-Function Block: '<S826>/B_82_302' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 1514, SS_CALL_MDL_OUTPUTS);

   
  
    if (ssGetTaskTime(S, 2) >= _rtB->B_97_1513_0) {
        _rtB->B_97_1521_0[0] = (_rtB->B_97_1512_0[0] - _rtB->B_97_1514_0[0]) * _rtP->P_1645 + (_rtP->P_496 * _rtB->B_97_1511_0[0] - _rtP->P_495 * _rtDW->UnitDelay_DSTATE_h[0]);
        _rtB->B_97_1521_0[1] = (_rtB->B_97_1512_0[1] - _rtB->B_97_1514_0[1]) * _rtP->P_1645 + (_rtP->P_496 * _rtB->B_97_1511_0[1] - _rtP->P_495 * _rtDW->UnitDelay_DSTATE_h[1]);
        _rtB->B_97_1521_0[2] = (_rtB->B_97_1512_0[2] - _rtB->B_97_1514_0[2]) * _rtP->P_1645 + (_rtP->P_496 * _rtB->B_97_1511_0[2] - _rtP->P_495 * _rtDW->UnitDelay_DSTATE_h[2]);
    } else {
        _rtB->B_97_1521_0[0] = _rtDW->UnitDelay1_DSTATE_k[0];
        _rtB->B_97_1521_0[1] = _rtDW->UnitDelay1_DSTATE_k[1];
        _rtB->B_97_1521_0[2] = _rtDW->UnitDelay1_DSTATE_k[2];
    }
    if (_rtDW->systemEnable_p != 0) {
        _rtDW->lastSin_f = muDoubleScalarSin(_rtP->P_1650 * ssGetTaskTime(S, 2));
        _rtDW->lastCos_g = muDoubleScalarCos(_rtP->P_1650 * ssGetTaskTime(S, 2));
        _rtDW->systemEnable_p = 0;
    }
    rtb_B_97_967_0 = ((_rtDW->lastSin_f * _rtP->P_1654 + _rtDW->lastCos_g * _rtP->P_1653) * _rtP->P_1652 + (_rtDW->lastCos_g * _rtP->P_1654 - _rtDW->lastSin_f * _rtP->P_1653) * _rtP->P_1651) * _rtP->P_1648 + _rtP->P_1649;
    _rtB->B_97_1523_0[0] = _rtB->B_97_265_0[0] * rtb_B_97_967_0;
    _rtB->B_97_1523_0[1] = _rtB->B_97_265_0[1] * rtb_B_97_967_0;
    _rtB->B_97_1523_0[2] = _rtB->B_97_265_0[2] * rtb_B_97_967_0;
    if (_rtDW->Integ4_SYSTEM_ENABLE_lp != 0) {
        _rtB->B_97_1524_0[0] = _rtDW->Integ4_DSTATE_hb[0];
        _rtB->B_97_1524_0[1] = _rtDW->Integ4_DSTATE_hb[1];
        _rtB->B_97_1524_0[2] = _rtDW->Integ4_DSTATE_hb[2];
    } else {
        _rtB->B_97_1524_0[0] = _rtP->P_1655 * _rtB->B_97_1523_0[0] + _rtDW->Integ4_DSTATE_hb[0];
        _rtB->B_97_1524_0[1] = _rtP->P_1655 * _rtB->B_97_1523_0[1] + _rtDW->Integ4_DSTATE_hb[1];
        _rtB->B_97_1524_0[2] = _rtP->P_1655 * _rtB->B_97_1523_0[2] + _rtDW->Integ4_DSTATE_hb[2];
    }
    _rtB->B_97_1525_0 = _rtP->P_1657;
                
          /* Level2 S-Function Block: '<S824>/B_82_304' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 1526, SS_CALL_MDL_OUTPUTS);

   
  
    if (ssGetTaskTime(S, 2) >= _rtB->B_97_1525_0) {
        _rtB->B_97_1533_0[0] = (_rtB->B_97_1524_0[0] - _rtB->B_97_1526_0[0]) * _rtP->P_1666 + (_rtP->P_494 * _rtB->B_97_1523_0[0] - _rtP->P_493 * _rtDW->UnitDelay_DSTATE_l3k[0]);
        _rtB->B_97_1533_0[1] = (_rtB->B_97_1524_0[1] - _rtB->B_97_1526_0[1]) * _rtP->P_1666 + (_rtP->P_494 * _rtB->B_97_1523_0[1] - _rtP->P_493 * _rtDW->UnitDelay_DSTATE_l3k[1]);
        _rtB->B_97_1533_0[2] = (_rtB->B_97_1524_0[2] - _rtB->B_97_1526_0[2]) * _rtP->P_1666 + (_rtP->P_494 * _rtB->B_97_1523_0[2] - _rtP->P_493 * _rtDW->UnitDelay_DSTATE_l3k[2]);
    } else {
        _rtB->B_97_1533_0[0] = _rtDW->UnitDelay1_DSTATE_d[0];
        _rtB->B_97_1533_0[1] = _rtDW->UnitDelay1_DSTATE_d[1];
        _rtB->B_97_1533_0[2] = _rtDW->UnitDelay1_DSTATE_d[2];
    }
    B_97_726_1_idx_0 = B_97_726_1_idx_0 * muDoubleScalarHypot(_rtB->B_97_1521_0[0], _rtB->B_97_1533_0[0]) * _rtP->P_1669;
    B_97_719_0_idx_0 = (_rtP->P_1670 * B_97_719_0_idx_0 - _rtP->P_1671 * muDoubleScalarAtan2(_rtB->B_97_1533_0[0], _rtB->B_97_1521_0[0])) * _rtP->P_1672;
    B_97_726_1_idx_1 = B_97_726_1_idx_1 * muDoubleScalarHypot(_rtB->B_97_1521_0[1], _rtB->B_97_1533_0[1]) * _rtP->P_1669;
    B_97_719_0_idx_1 = (_rtP->P_1670 * B_97_719_0_idx_1 - _rtP->P_1671 * muDoubleScalarAtan2(_rtB->B_97_1533_0[1], _rtB->B_97_1521_0[1])) * _rtP->P_1672;
    B_97_726_1_idx_2 = B_97_726_1_idx_2 * muDoubleScalarHypot(_rtB->B_97_1521_0[2], _rtB->B_97_1533_0[2]) * _rtP->P_1669;
    rtb_B_97_1_0 = (_rtP->P_1670 * B_97_719_0_idx_2 - _rtP->P_1671 * muDoubleScalarAtan2(_rtB->B_97_1533_0[2], _rtB->B_97_1521_0[2])) * _rtP->P_1672;
    muDoubleScalarSinCos(B_97_719_0_idx_0, &B_97_1542_0_idx_0, &rtb_B_97_1173_0);
    muDoubleScalarSinCos(B_97_719_0_idx_1, &B_97_1542_0_idx_1, &B_97_1542_1_idx_1);
    muDoubleScalarSinCos(rtb_B_97_1_0, &B_97_1542_0_idx_2, &B_97_1542_1_idx_2);
    _rtB->B_97_1544_0 = (B_97_726_1_idx_0 * rtb_B_97_1173_0 + B_97_726_1_idx_1 * B_97_1542_1_idx_1) + B_97_726_1_idx_2 * B_97_1542_1_idx_2;
}
_rtB->B_97_1545_0 = _rtX->integ1_CSTATE_do;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T_PWORK_dy.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T_PWORK_dy.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_1674;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        B_97_1546_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T_IWORK_i.CircularBufSize,
        &_rtDW->T_IWORK_i.Last,
        _rtDW->T_IWORK_i.Tail,
        _rtDW->T_IWORK_i.Head,
        _rtP->P_1675,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
B_97_1551_0.re = _rtB->B_97_1545_0 - B_97_1546_0;
_rtB->B_97_1548_0 = _rtX->Integ2_CSTATE_oo;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK_ck.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK_ck.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_1677;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        B_97_1549_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T1_IWORK_g.CircularBufSize,
        &_rtDW->T1_IWORK_g.Last,
        _rtDW->T1_IWORK_g.Tail,
        _rtDW->T1_IWORK_g.Head,
        _rtP->P_1678,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
B_97_1551_0.im = _rtB->B_97_1548_0 - B_97_1549_0;
_rtB->B_97_1552_0 = _rtX->integ1_CSTATE_m;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T_PWORK_p.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T_PWORK_p.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_1680;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        B_97_1553_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T_IWORK_i5.CircularBufSize,
        &_rtDW->T_IWORK_i5.Last,
        _rtDW->T_IWORK_i5.Tail,
        _rtDW->T_IWORK_i5.Head,
        _rtP->P_1681,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
B_97_1558_0.re = _rtB->B_97_1552_0 - B_97_1553_0;
_rtB->B_97_1555_0 = _rtX->Integ2_CSTATE_ot;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK_ph.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK_ph.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_1683;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        B_97_1556_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T1_IWORK_i.CircularBufSize,
        &_rtDW->T1_IWORK_i.Last,
        _rtDW->T1_IWORK_i.Tail,
        _rtDW->T1_IWORK_i.Head,
        _rtP->P_1684,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
B_97_1558_0.im = _rtB->B_97_1555_0 - B_97_1556_0;
_rtB->B_97_1559_0 = _rtX->integ1_CSTATE_i;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T_PWORK_bl.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T_PWORK_bl.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_1686;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        B_97_1560_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T_IWORK_by.CircularBufSize,
        &_rtDW->T_IWORK_by.Last,
        _rtDW->T_IWORK_by.Tail,
        _rtDW->T_IWORK_by.Head,
        _rtP->P_1687,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
B_97_1565_0.re = _rtB->B_97_1559_0 - B_97_1560_0;
_rtB->B_97_1562_0 = _rtX->Integ2_CSTATE_j;
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK_oh.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK_oh.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_1689;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        B_97_1563_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->T1_IWORK_p.CircularBufSize,
        &_rtDW->T1_IWORK_p.Last,
        _rtDW->T1_IWORK_p.Tail,
        _rtDW->T1_IWORK_p.Head,
        _rtP->P_1690,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
B_97_1565_0.im = _rtB->B_97_1562_0 - B_97_1563_0;
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    _rtB->B_97_1567_0 = _rtB->B_97_1566_0;
}





      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_NegSeqComputation(S, _rtB->B_97_1567_0, B_97_1551_0, B_97_1558_0, B_97_1565_0, &_rtB->PosSeqComputation_g, &_rtDW->PosSeqComputation_g, &_rtP->PosSeqComputation_g);
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    _rtB->B_97_1570_0 = _rtB->B_97_1569_0;
}





      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_NegSeqComputation(S, _rtB->B_97_1570_0, B_97_1551_0, B_97_1558_0, B_97_1565_0, &_rtB->NegSeqComputation_e, &_rtDW->NegSeqComputation_e, &_rtP->NegSeqComputation_e);
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    _rtB->B_97_1573_0 = _rtB->B_97_1572_0;
}





      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_ZeroSeqComputation(S, _rtB->B_97_1573_0, B_97_1551_0, B_97_1558_0, B_97_1565_0, &_rtB->ZeroSeqComputation_k, &_rtDW->ZeroSeqComputation_k, &_rtP->ZeroSeqComputation_k);
_rtB->B_97_1575_0[0] = muDoubleScalarHypot(_rtB->PosSeqComputation_g.B_39_2_0.re, _rtB->PosSeqComputation_g.B_39_2_0.im);
_rtB->B_97_1575_0[1] = muDoubleScalarHypot(_rtB->NegSeqComputation_e.B_39_2_0.re, _rtB->NegSeqComputation_e.B_39_2_0.im);
_rtB->B_97_1575_0[2] = muDoubleScalarHypot(_rtB->ZeroSeqComputation_k.B_41_1_0.re, _rtB->ZeroSeqComputation_k.B_41_1_0.im);
rtb_B_97_30_0[0] = muDoubleScalarAtan2(_rtB->PosSeqComputation_g.B_39_2_0.im, _rtB->PosSeqComputation_g.B_39_2_0.re);
rtb_B_97_30_0[1] = muDoubleScalarAtan2(_rtB->NegSeqComputation_e.B_39_2_0.im, _rtB->NegSeqComputation_e.B_39_2_0.re);
rtb_B_97_30_0[2] = muDoubleScalarAtan2(_rtB->ZeroSeqComputation_k.B_41_1_0.im, _rtB->ZeroSeqComputation_k.B_41_1_0.re);
_rtB->B_97_1576_0 = _rtP->P_1694 * _rtB->B_97_1575_0[0];
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_1577_0 = _rtDW->UnitDelay_DSTATE_jz;
    B_97_726_1_idx_0 *= B_97_1542_0_idx_0;
    B_97_726_1_idx_1 *= B_97_1542_0_idx_1;
    _rtB->B_97_1579_0 = (B_97_726_1_idx_0 + B_97_726_1_idx_1) + B_97_726_1_idx_2 * B_97_1542_0_idx_2;
}
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 1580, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 1581, SS_CALL_MDL_OUTPUTS);
   
  
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 1582, SS_CALL_MDL_OUTPUTS);
   
  
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 1583, SS_CALL_MDL_OUTPUTS);
   
  
}
_rtB->B_97_1585_0 = _rtP->P_1697 * _rtB->B_97_1575_0[0];
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    if (ssIsMajorTimeStep(S) != 0) {
        _rtDW->RelationalOperator_Mode_c[0] = (_rtB->B_97_1575_0[1] > _rtB->B_97_1585_0);
        _rtDW->RelationalOperator_Mode_c[1] = (_rtB->B_97_1575_0[2] > _rtB->B_97_1585_0);
    }
    _rtB->B_97_1586_0[0] = _rtDW->RelationalOperator_Mode_c[0];
    _rtB->B_97_1586_0[1] = _rtDW->RelationalOperator_Mode_c[1];
}
if (_rtB->B_97_1586_0[0] >= _rtP->P_1698) {
    B_97_1587_0[0] = rtb_B_97_30_0[1];
} else {
    B_97_1587_0[0] = _rtB->B_97_1584_0;
}
if (_rtB->B_97_1586_0[1] >= _rtP->P_1698) {
    B_97_1587_0[1] = rtb_B_97_30_0[2];
} else {
    B_97_1587_0[1] = _rtB->B_97_1584_0;
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_1588_0 = _rtP->P_1699 * _rtB->B_97_295_0[0];
}
rtb_B_97_967_0 = muDoubleScalarSin(_rtP->P_1702 * ssGetTaskTime(S, 0) + _rtP->P_1703) * _rtP->P_1700 + _rtP->P_1701;
_rtB->B_97_1590_0 = rtb_B_97_967_0 * _rtB->B_97_1588_0;
B_97_1331_0_idx_0 = muDoubleScalarSin(_rtP->P_1706 * ssGetTaskTime(S, 0) + _rtP->P_1707) * _rtP->P_1704 + _rtP->P_1705;
_rtB->B_97_1592_0 = _rtB->B_97_1588_0 * B_97_1331_0_idx_0;
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_1593_0 = _rtP->P_1708 * _rtB->B_97_295_0[1];
}
_rtB->B_97_1594_0 = rtb_B_97_967_0 * _rtB->B_97_1593_0;
_rtB->B_97_1595_0 = _rtB->B_97_1593_0 * B_97_1331_0_idx_0;
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_1596_0 = _rtP->P_1709 * _rtB->B_97_295_0[2];
}
_rtB->B_97_1597_0 = rtb_B_97_967_0 * _rtB->B_97_1596_0;
_rtB->B_97_1598_0 = _rtB->B_97_1596_0 * B_97_1331_0_idx_0;
B_97_1599_0 = _rtP->P_1710 * rtb_B_97_30_0[0];
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    B_97_1542_0_idx_0 = _rtP->P_1711 * _rtB->B_97_295_0[0];
    B_97_1542_0_idx_1 = _rtP->P_1711 * _rtB->B_97_295_0[1];
    B_97_1542_0_idx_2 = _rtP->P_1711 * _rtB->B_97_295_0[2];
    B_97_1331_0_idx_0 = muDoubleScalarMod(_rtDW->DiscreteTimeIntegrator_DSTATE_om, _rtB->B_97_1602_0);
    rtb_B_97_967_0 = muDoubleScalarSin(B_97_1331_0_idx_0);
    B_97_1331_0_idx_0 = muDoubleScalarCos(B_97_1331_0_idx_0);
    rtb_B_97_16_0 = (0.0 - rtb_B_97_967_0 * _rtB->B_97_1606_0) - B_97_1331_0_idx_0 * _rtB->B_97_1608_0;
    rtb_B_97_17_0 = rtb_B_97_967_0 * _rtB->B_97_1608_0 - B_97_1331_0_idx_0 * _rtB->B_97_1606_0;
    rtb_B_97_1236_0 = (0.0 - rtb_B_97_16_0) - rtb_B_97_967_0;
    rtb_B_97_1_0 = (0.0 - rtb_B_97_17_0) - B_97_1331_0_idx_0;
    _rtB->B_97_1616_0 = ((B_97_1542_0_idx_0 * B_97_1331_0_idx_0 + B_97_1542_0_idx_1 * rtb_B_97_17_0) + B_97_1542_0_idx_2 * rtb_B_97_1_0) * 0.66666666666666663;
    if (_rtDW->Integ4_SYSTEM_ENABLE_b != 0) {
        _rtB->B_97_1617_0 = _rtDW->Integ4_DSTATE_m2;
    } else {
        _rtB->B_97_1617_0 = _rtP->P_1717 * _rtB->B_97_1616_0 + _rtDW->Integ4_DSTATE_m2;
    }
    if (_rtB->B_97_1577_0 > _rtP->P_1719) {
        B_5_0_0_idx_3 = _rtP->P_1719;
    } else if (_rtB->B_97_1577_0 < _rtP->P_1720) {
        B_5_0_0_idx_3 = _rtP->P_1720;
    } else {
        B_5_0_0_idx_3 = _rtB->B_97_1577_0;
    }
    rtb_B_97_1173_0 = 1.0 / B_5_0_0_idx_3 / 1.0e-5;
    rtb_B_97_5_0 = muDoubleScalarCeil(rtb_B_97_1173_0);
    _rtB->B_97_1621_0 = _rtP->P_1721 * rtb_B_97_5_0;
                
          /* Level2 S-Function Block: '<S812>/B_82_319' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 1622, SS_CALL_MDL_OUTPUTS);

   
  
    rtb_B_97_82_0 = ssGetTaskTime(S, 2);
    if (rtb_B_97_82_0 < _rtP->P_1723) {
        rtb_B_97_82_0 = _rtP->P_1724;
    } else {
        rtb_B_97_82_0 = _rtP->P_1725;
    }
    if (rtb_B_97_82_0 >= _rtP->P_1727) {
        rtb_B_97_5_0 = rtb_B_97_1173_0 - rtb_B_97_5_0;
        rtb_B_97_1173_0 = ((_rtB->B_97_1616_0 - _rtDW->UnitDelay_DSTATE_o0) * rtb_B_97_5_0 * _rtP->P_488 + _rtB->B_97_1616_0) * (rtb_B_97_5_0 / rtb_B_97_1173_0) + (_rtB->B_97_1617_0 - _rtB->B_97_1622_0) * _rtB->B_97_1577_0;
    } else {
        rtb_B_97_1173_0 = _rtP->P_1726;
    }
    _rtB->B_97_1629_0 = _rtB->B_97_1628_0;
    i = ssIsSampleHit(S, 2, 0);
    if ((i != 0) && (ssIsMajorTimeStep(S) != 0)) {
        if (_rtB->B_97_1629_0 > 0.0) {
            if (!_rtDW->AutomaticGainControl_MODE_o) {
                if (ssGetTaskTime(S, 2) != ssGetTStart(S)) {
                    ssSetBlockStateForSolverChangedAtMajorStep(S);
                }
  
                _rtDW->Integ4_SYSTEM_ENABLE_au = 1U;
                _rtDW->Integ4_SYSTEM_ENABLE_d = 1U;
                _rtDW->AutomaticGainControl_MODE_o = true;
            }
        } else {
            if (_rtDW->AutomaticGainControl_MODE_o) {
                ssSetBlockStateForSolverChangedAtMajorStep(S);
  
                _rtDW->AutomaticGainControl_MODE_o = false;
            }
        }
    }
    if (_rtDW->AutomaticGainControl_MODE_o) {
        i = ssIsSampleHit(S, 2, 0);
        if (i != 0) {
            rtb_B_97_5_0 = (0.0 - rtb_B_97_967_0 * _rtB->B_77_0_0) - B_97_1331_0_idx_0 * _rtB->B_77_2_0;
            rtb_B_97_5_1 = rtb_B_97_967_0 * _rtB->B_77_2_0 - B_97_1331_0_idx_0 * _rtB->B_77_0_0;
            rtb_B_97_27_0 = (0.0 - rtb_B_97_5_0) - rtb_B_97_967_0;
            rtb_B_97_28_0 = (0.0 - rtb_B_97_5_1) - B_97_1331_0_idx_0;
            _rtB->B_77_10_0 = ((B_97_1542_0_idx_0 * rtb_B_97_967_0 + B_97_1542_0_idx_1 * rtb_B_97_5_0) + B_97_1542_0_idx_2 * rtb_B_97_27_0) * 0.66666666666666663;
            if (_rtDW->Integ4_SYSTEM_ENABLE_au != 0) {
                _rtB->B_77_11_0 = _rtDW->Integ4_DSTATE_fi;
            } else {
                _rtB->B_77_11_0 = _rtP->P_462 * _rtB->B_77_10_0 + _rtDW->Integ4_DSTATE_fi;
            }
            if (_rtB->B_97_1577_0 > _rtP->P_464) {
                B_5_0_0_idx_3 = _rtP->P_464;
            } else if (_rtB->B_97_1577_0 < _rtP->P_465) {
                B_5_0_0_idx_3 = _rtP->P_465;
            } else {
                B_5_0_0_idx_3 = _rtB->B_97_1577_0;
            }
            B_97_726_1_idx_0 = 1.0 / B_5_0_0_idx_3 / 1.0e-5;
            B_97_726_1_idx_1 = muDoubleScalarCeil(B_97_726_1_idx_0);
            _rtB->B_77_15_0 = _rtP->P_466 * B_97_726_1_idx_1;
                
          /* Level2 S-Function Block: '<S807>/B_66_0' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 77, 16, SS_CALL_MDL_OUTPUTS);

   
  
            rtb_B_97_82_0 = ssGetTaskTime(S, 2);
            if (rtb_B_97_82_0 < _rtP->P_468) {
                rtb_B_97_82_0 = _rtP->P_469;
            } else {
                rtb_B_97_82_0 = _rtP->P_470;
            }
            if (rtb_B_97_82_0 >= _rtP->P_472) {
                B_97_726_1_idx_1 = B_97_726_1_idx_0 - B_97_726_1_idx_1;
                B_97_726_1_idx_0 = ((_rtB->B_77_10_0 - _rtDW->UnitDelay_DSTATE_hk) * B_97_726_1_idx_1 * _rtP->P_457 + _rtB->B_77_10_0) * (B_97_726_1_idx_1 / B_97_726_1_idx_0) + (_rtB->B_77_11_0 - _rtB->B_77_16_0) * _rtB->B_97_1577_0;
            } else {
                B_97_726_1_idx_0 = _rtP->P_471;
            }
            _rtB->B_77_22_0 = ((B_97_1542_0_idx_0 * B_97_1331_0_idx_0 + B_97_1542_0_idx_1 * rtb_B_97_5_1) + B_97_1542_0_idx_2 * rtb_B_97_28_0) * 0.66666666666666663;
            if (_rtDW->Integ4_SYSTEM_ENABLE_d != 0) {
                _rtB->B_77_23_0 = _rtDW->Integ4_DSTATE_ez;
            } else {
                _rtB->B_77_23_0 = _rtP->P_473 * _rtB->B_77_22_0 + _rtDW->Integ4_DSTATE_ez;
            }
            if (_rtB->B_97_1577_0 > _rtP->P_475) {
                B_5_0_0_idx_3 = _rtP->P_475;
            } else if (_rtB->B_97_1577_0 < _rtP->P_476) {
                B_5_0_0_idx_3 = _rtP->P_476;
            } else {
                B_5_0_0_idx_3 = _rtB->B_97_1577_0;
            }
            B_97_726_1_idx_1 = 1.0 / B_5_0_0_idx_3 / 1.0e-5;
            rtb_B_97_966_0 = muDoubleScalarCeil(B_97_726_1_idx_1);
            _rtB->B_77_27_0 = _rtP->P_477 * rtb_B_97_966_0;
                
          /* Level2 S-Function Block: '<S809>/B_66_2' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 77, 28, SS_CALL_MDL_OUTPUTS);

   
  
            rtb_B_97_82_0 = ssGetTaskTime(S, 2);
            if (rtb_B_97_82_0 < _rtP->P_479) {
                rtb_B_97_82_0 = _rtP->P_480;
            } else {
                rtb_B_97_82_0 = _rtP->P_481;
            }
            if (rtb_B_97_82_0 >= _rtP->P_483) {
                rtb_B_97_966_0 = B_97_726_1_idx_1 - rtb_B_97_966_0;
                B_97_726_1_idx_1 = ((_rtB->B_77_22_0 - _rtDW->UnitDelay_DSTATE_fvk) * rtb_B_97_966_0 * _rtP->P_458 + _rtB->B_77_22_0) * (rtb_B_97_966_0 / B_97_726_1_idx_1) + (_rtB->B_77_23_0 - _rtB->B_77_28_0) * _rtB->B_97_1577_0;
            } else {
                B_97_726_1_idx_1 = _rtP->P_482;
            }
            B_77_36_0 = _rtP->P_484 * muDoubleScalarAtan2(B_97_726_1_idx_1, B_97_726_1_idx_0);
            tmpForInput[0] = B_97_1542_0_idx_0;
            tmpForInput[1] = B_97_1542_0_idx_1;
            tmpForInput[2] = B_97_1542_0_idx_2;
            tmpForInput[3] = rtb_B_97_967_0;
            tmpForInput[4] = B_97_1331_0_idx_0;
            tmpForInput[5] = rtb_B_97_5_0;
            tmpForInput[6] = rtb_B_97_5_1;
            tmpForInput[7] = rtb_B_97_27_0;
            tmpForInput[8] = rtb_B_97_28_0;
            rtb_B_77_37_0 = -0.0;
            for (i = 0; i < 9; i++) {
                rtb_B_77_37_0 += tmpForInput[i];
            }
            B_77_38_0 = _rtP->P_485 * rtb_B_77_37_0;
            B_5_0_0_idx_3 = muDoubleScalarHypot(B_97_726_1_idx_0, B_97_726_1_idx_1);
            if (B_5_0_0_idx_3 > _rtP->P_486) {
                B_5_0_0_idx_3 = _rtP->P_486;
            } else {
                if (B_5_0_0_idx_3 < _rtP->P_487) {
                    B_5_0_0_idx_3 = _rtP->P_487;
                }
            }
            _rtB->B_77_40_0 = 1.0 / B_5_0_0_idx_3;
        }
        if (ssIsMajorTimeStep(S) != 0) {
            srUpdateBC(_rtDW->AutomaticGainControl_SubsysRanBC_h);
        }
    }
    rtb_B_77_37_0 = rtb_B_97_1173_0 * _rtB->B_77_40_0;
    _rtB->B_97_1636_0 = _rtP->P_1734 * rtb_B_77_37_0 * _rtP->P_1735;
    B_97_1637_0 = _rtDW->UD_DSTATE_k;
    B_97_1638_0 = _rtB->B_97_1636_0 - B_97_1637_0;
    B_5_0_0_idx_3 = (_rtP->P_1729 * rtb_B_77_37_0 + _rtDW->DiscreteTimeIntegrator_DSTATE_k) + B_97_1638_0;
    if (B_5_0_0_idx_3 > _rtP->P_1737) {
        _rtB->B_97_1640_0 = _rtP->P_1737;
    } else if (B_5_0_0_idx_3 < _rtP->P_1738) {
        _rtB->B_97_1640_0 = _rtP->P_1738;
    } else {
        _rtB->B_97_1640_0 = B_5_0_0_idx_3;
    }
    B_5_0_0_idx_3 = _rtP->P_1739 * _rtB->B_97_1640_0 - _rtDW->UnitDelay_DSTATE_d0;
    if (B_5_0_0_idx_3 > _rtP->P_1741) {
        B_5_0_0_idx_3 = _rtP->P_1741;
    } else {
        if (B_5_0_0_idx_3 < _rtP->P_1742) {
            B_5_0_0_idx_3 = _rtP->P_1742;
        }
    }
    _rtB->B_97_1645_0 = B_5_0_0_idx_3 + _rtDW->UnitDelay_DSTATE_d0;
          {
                 _rtB->B_97_1646_0 = (_rtP->P_1745[0])*_rtDW->DiscreteStateSpace_DSTATE_e[0] + (_rtP->P_1745[1])*_rtDW->DiscreteStateSpace_DSTATE_e[1];
    _rtB->B_97_1646_0 += _rtP->P_1746*_rtB->B_97_1645_0;
  }

   
  
    _rtB->B_97_1648_0 = _rtP->P_1748 * rtb_B_77_37_0;
    B_97_1649_0 = ((B_97_1542_0_idx_0 * rtb_B_97_967_0 + B_97_1542_0_idx_1 * rtb_B_97_16_0) + B_97_1542_0_idx_2 * rtb_B_97_1236_0) * 0.66666666666666663;
    tmpForInput[0] = B_97_1542_0_idx_0;
    tmpForInput[1] = B_97_1542_0_idx_1;
    tmpForInput[2] = B_97_1542_0_idx_2;
    tmpForInput[3] = rtb_B_97_967_0;
    tmpForInput[4] = B_97_1331_0_idx_0;
    tmpForInput[5] = rtb_B_97_16_0;
    tmpForInput[6] = rtb_B_97_17_0;
    tmpForInput[7] = rtb_B_97_1236_0;
    tmpForInput[8] = rtb_B_97_1_0;
    B_97_1542_0_idx_0 = -0.0;
    for (i = 0; i < 9; i++) {
        B_97_1542_0_idx_0 += tmpForInput[i];
    }
    B_97_1651_0 = _rtP->P_1749 * B_97_1542_0_idx_0;
}
          {


      
    
    


        
        real_T **uBuffer = (real_T**)&_rtDW->ENGINETd_PWORK_d.TUbufferPtrs[0];
        real_T **tBuffer = (real_T**)&_rtDW->ENGINETd_PWORK_d.TUbufferPtrs[1];
        real_T simTime   = ssGetT(S);   

        real_T tMinusDelay = simTime - _rtP->P_1750;    
    
       
    
    
  
    
      
                     


          
        

      
      
        
        B_97_1714_0 = Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        _rtDW->ENGINETd_IWORK_c.CircularBufSize,
        &_rtDW->ENGINETd_IWORK_c.Last,
        _rtDW->ENGINETd_IWORK_c.Tail,
        _rtDW->ENGINETd_IWORK_c.Head,
        _rtP->P_1751,
        0,
        (boolean_T) (ssIsMinorTimeStep(S) && (ssGetTimeOfLastOutput(S) == ssGetT(S))));  	    
      
      
  }
   
  
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    rtb_B_97_1177_0[0] = _rtP->P_1753[0] * rtb_B_97_56_0[0] * _rtB->B_97_55_0[1];
    rtb_B_97_1177_0[1] = _rtP->P_1753[1] * rtb_B_97_56_0[1] * _rtB->B_97_55_0[0];
    rtb_B_97_1718_0 = rtb_B_97_1177_0[0] + rtb_B_97_1177_0[1];
    i = ssIsSampleHit(S, 2, 0);
    if (i != 0) {
        _rtB->B_97_1725_0[0] = rtb_B_97_32_0;
        _rtB->B_97_1725_0[1] = _rtP->P_1752 * rtb_B_97_48_0;
        _rtB->B_97_1725_0[2] = rtb_B_97_48_0 * rtb_B_97_1718_0 * _rtP->P_1754;
        _rtB->B_97_1725_0[3] = rtb_B_97_47_0;
        _rtB->B_97_1725_0[4] = _rtP->P_1757 * muDoubleScalarRem(_rtDW->DiscreteTimeIntegrator2_DSTATE_j, 6.2831853071795862);
        _rtB->B_97_1725_0[5] = _rtP->P_1758 * rtb_B_97_1718_0;
    }
}
_rtB->B_97_1726_0 = B_97_1714_0 * _rtB->B_97_1725_0[1];
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_1727_0 = _rtDW->UnitDelay2_DSTATE_h;
    i = ssIsSampleHit(S, 2, 0);
    if (i != 0) {
        rtb_B_97_1740_0[0] = _rtB->B_97_116_0[27];
        rtb_B_97_1740_0[1] = _rtB->B_97_116_0[28];
    }
    rtb_B_97_1740_0[0] *= _rtP->P_1766;
    rtb_B_97_1740_0[1] *= _rtP->P_1766;
    _rtB->B_97_1742_0 = ((2.0 * rtb_B_97_1740_0[0] + rtb_B_97_1740_0[1]) * rtb_B_97_36_1 + 1.7320508075688772 * rtb_B_97_1740_0[1] * rtb_B_97_36_0) * 0.33333333333333331;
    _rtB->B_97_1743_0 = ((2.0 * rtb_B_97_1740_0[0] + rtb_B_97_1740_0[1]) * rtb_B_97_36_0 + -1.7320508075688772 * rtb_B_97_1740_0[1] * rtb_B_97_36_1) * 0.33333333333333331;
    rtb_B_97_1194_0[11] = _rtP->P_1773[11] * _rtB->B_97_1742_0;
    rtb_B_97_1194_0[12] = _rtP->P_1773[12] * _rtB->B_97_1743_0;
    rtb_B_97_67_0 = muDoubleScalarPower(rtb_B_97_1194_0[12], 2.0) + muDoubleScalarPower(rtb_B_97_1194_0[11], 2.0);
    if (rtb_B_97_67_0 < 0.0) {
        _rtB->B_97_1759_0 = -muDoubleScalarSqrt(-rtb_B_97_67_0);
    } else {
        _rtB->B_97_1759_0 = muDoubleScalarSqrt(rtb_B_97_67_0);
    }
}
            /* Call into Simulink for Scope */
    ssCallAccelRunBlock(S, 97, 1760, SS_CALL_MDL_OUTPUTS);
   
  
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    rtb_B_97_67_0 = muDoubleScalarPower(rtb_B_97_1194_0[12], 2.0) + muDoubleScalarPower(rtb_B_97_1194_0[11], 2.0);
    if (rtb_B_97_67_0 < 0.0) {
        _rtB->B_97_1763_0 = -muDoubleScalarSqrt(-rtb_B_97_67_0);
    } else {
        _rtB->B_97_1763_0 = muDoubleScalarSqrt(rtb_B_97_67_0);
    }
          {
                 B_97_1764_0 = _rtP->P_1778*_rtDW->DiscreteStateSpace_DSTATE_ax;
    B_97_1764_0 += _rtP->P_1779*_rtB->B_97_1763_0;
  }

   
  
    _rtB->B_97_1766_0 = ((_rtB->B_97_1761_0 + _rtB->B_97_1762_0) - B_97_1764_0) - _rtDW->UnitDelay1_DSTATE_fx;
          {
                 B_97_1767_0 = _rtP->P_1785*_rtB->B_97_1766_0;
  }

   
  
    _rtB->B_97_1768_0 = _rtP->P_1787 * B_97_1767_0;
          {
                 B_97_1769_0 = _rtP->P_1790*_rtDW->DiscreteStateSpace_DSTATE_m5;
    B_97_1769_0 += _rtP->P_1791*_rtB->B_97_1768_0;
  }

   
  
    rtb_B_97_36_0 = (B_97_1769_0 > _rtB->B_97_1770_0);
    rtb_B_97_36_1 = _rtP->P_1794 * B_97_1764_0 + _rtB->B_97_1777_0;
    rtb_B_97_32_0 = (B_97_1769_0 < rtb_B_97_36_1);
    _rtB->B_97_1790_0 = (((real_T)((rtb_B_97_36_0 != 0.0) && (rtb_B_97_32_0 != 0.0)) * B_97_1769_0 + (real_T)!(rtb_B_97_36_0 != 0.0) * _rtB->B_97_1770_0) + (real_T)!(rtb_B_97_32_0 != 0.0) * rtb_B_97_36_1) * _rtP->P_1798;
          {
                 _rtB->B_97_1791_0 = _rtP->P_1801*_rtDW->DiscreteStateSpace_DSTATE_l;
    _rtB->B_97_1791_0 += _rtP->P_1802*_rtB->B_97_1790_0;
  }

   
  
    _rtB->B_97_1792_0 = _rtP->P_1804 * _rtB->B_97_1791_0;
          {
                 _rtB->B_97_1793_0 = (_rtP->P_1807)*_rtDW->DiscreteStateSpace_DSTATE_hk;
    _rtB->B_97_1793_0 += _rtP->P_1808*_rtB->B_97_1792_0;
  }

   
  
}
_rtB->B_97_1794_0 = 0.0;
_rtB->B_97_1794_0 += _rtP->P_1811[0] * _rtX->CONTROLSYSTEM_CSTATE_b[0];
_rtB->B_97_1794_0 += _rtP->P_1811[1] * _rtX->CONTROLSYSTEM_CSTATE_b[1];
_rtB->B_97_1795_0 = _rtP->P_1812 * _rtB->B_97_1794_0;
/* Limited  Integrator  */
if (ssIsMajorTimeStep(S) != 0) {
    if (_rtX->Integrator_CSTATE_k >= _rtP->P_1814) {
        _rtX->Integrator_CSTATE_k = _rtP->P_1814;
    } else {
        if (_rtX->Integrator_CSTATE_k <= _rtP->P_1815) {
            _rtX->Integrator_CSTATE_k = _rtP->P_1815;
        }
    }
}
_rtB->B_97_1796_0 = _rtX->Integrator_CSTATE_k;
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtB->B_97_1798_0 = _rtB->B_97_1797_0 - _rtB->B_97_1725_0[1];
}
_rtB->B_97_1799_0 = 0.0;
_rtB->B_97_1799_0 += _rtP->P_1818 * _rtX->TF1_CSTATE_p;
_rtB->B_97_1799_0 += _rtP->P_1819 * _rtB->B_97_1795_0;
_rtB->B_97_1800_0 = 0.0;
_rtB->B_97_1800_0 += _rtP->P_1821 * _rtX->TF2_CSTATE_m;
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    B_97_1801_0 = muDoubleScalarHypot(_rtB->B_97_1742_0, _rtB->B_97_1743_0);
    _rtB->B_97_1802_0 = _rtP->P_1822 * _rtB->B_97_1727_0;
}
i = ssIsSampleHit(S, 1, 0);
if ((i != 0) && (_rtDW->RateTransition1_semaphoreTaken_g == 0)) {
    _rtDW->RateTransition1_Buffer0_j = _rtB->B_97_1726_0;
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtDW->RateTransition1_semaphoreTaken_g = 1;
    rtb_B_97_1804_0 = _rtDW->RateTransition1_Buffer0_j;
    _rtDW->RateTransition1_semaphoreTaken_g = 0;
    rtb_B_97_1804_0 *= _rtP->P_1824;
    _rtB->B_97_1809_0 = ((rtb_B_97_1804_0 / (rtb_B_97_48_0 + 2.2204460492503131e-16) - rtb_B_97_1718_0) - _rtP->P_1825 * rtb_B_97_48_0) * _rtP->P_1826;
    _rtB->B_97_1810_0 = _rtP->P_1827 * rtb_B_97_47_0;
    _rtB->B_97_1811_0 = _rtP->P_1828 * rtb_B_97_48_0;
    if (_rtB->B_97_1896_0 >= _rtP->P_1834) {
        rtb_B_97_1804_0 = _rtB->B_97_1895_0;
    } else {
        rtb_B_97_1804_0 = rt_Lookup(_rtP->P_1832, 4, ssGetTaskTime(S, 2), _rtP->P_1833);
    }
    if (_rtB->B_97_1900_0 >= _rtP->P_1837) {
        _rtB->B_97_1903_0 = rtb_B_97_1804_0;
    } else {
        _rtB->B_97_1903_0 = _rtB->B_97_1901_0;
    }
}
_rtB->B_97_1905_0 = rt_Lookup(_rtP->P_1838, 5, ssGetT(S), _rtP->P_1839);
if (_rtB->B_97_1894_0 >= _rtP->P_1840) {
    _rtB->B_97_1906_0 = _rtB->B_97_1903_0;
} else {
    _rtB->B_97_1906_0 = _rtB->B_97_1905_0;
}
if ((_rtDW->TimeStampA >= ssGetT(S)) && (_rtDW->TimeStampB >= ssGetT(S))) {
    _rtB->B_97_1907_0 = 0.0;
} else {
    rtb_B_97_47_0 = _rtDW->TimeStampA;
    lastU = &_rtDW->LastUAtTimeA;
    if (_rtDW->TimeStampA < _rtDW->TimeStampB) {
        if (_rtDW->TimeStampB < ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB;
            lastU = &_rtDW->LastUAtTimeB;
        }
    } else {
        if (_rtDW->TimeStampA >= ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB;
            lastU = &_rtDW->LastUAtTimeB;
        }
    }
    _rtB->B_97_1907_0 = (_rtB->B_97_1905_0 - *lastU) / (ssGetT(S) - rtb_B_97_47_0);
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    if (_rtB->B_97_1913_0 >= _rtP->P_1844) {
        _rtB->B_97_1915_0 = rtb_B_97_1804_0;
    } else {
        _rtB->B_97_1915_0 = _rtB->B_97_1901_0;
    }
}
_rtB->B_97_1917_0 = rt_Lookup(_rtP->P_1845, 5, ssGetT(S), _rtP->P_1846);
if (_rtB->B_97_1912_0 >= _rtP->P_1847) {
    _rtB->B_97_1918_0 = _rtB->B_97_1915_0;
} else {
    _rtB->B_97_1918_0 = _rtB->B_97_1917_0;
}
if ((_rtDW->TimeStampA_c >= ssGetT(S)) && (_rtDW->TimeStampB_a >= ssGetT(S))) {
    _rtB->B_97_1919_0 = 0.0;
} else {
    rtb_B_97_47_0 = _rtDW->TimeStampA_c;
    lastU = &_rtDW->LastUAtTimeA_d;
    if (_rtDW->TimeStampA_c < _rtDW->TimeStampB_a) {
        if (_rtDW->TimeStampB_a < ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB_a;
            lastU = &_rtDW->LastUAtTimeB_d;
        }
    } else {
        if (_rtDW->TimeStampA_c >= ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB_a;
            lastU = &_rtDW->LastUAtTimeB_d;
        }
    }
    _rtB->B_97_1919_0 = (_rtB->B_97_1917_0 - *lastU) / (ssGetT(S) - rtb_B_97_47_0);
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    if (_rtB->B_97_1925_0 >= _rtP->P_1851) {
        _rtB->B_97_1927_0 = rtb_B_97_1804_0;
    } else {
        _rtB->B_97_1927_0 = _rtB->B_97_1901_0;
    }
}
_rtB->B_97_1929_0 = rt_Lookup(_rtP->P_1852, 5, ssGetT(S), _rtP->P_1853);
if (_rtB->B_97_1924_0 >= _rtP->P_1854) {
    _rtB->B_97_1930_0 = _rtB->B_97_1927_0;
} else {
    _rtB->B_97_1930_0 = _rtB->B_97_1929_0;
}
if ((_rtDW->TimeStampA_cb >= ssGetT(S)) && (_rtDW->TimeStampB_ax >= ssGetT(S))) {
    _rtB->B_97_1931_0 = 0.0;
} else {
    rtb_B_97_47_0 = _rtDW->TimeStampA_cb;
    lastU = &_rtDW->LastUAtTimeA_d2;
    if (_rtDW->TimeStampA_cb < _rtDW->TimeStampB_ax) {
        if (_rtDW->TimeStampB_ax < ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB_ax;
            lastU = &_rtDW->LastUAtTimeB_m;
        }
    } else {
        if (_rtDW->TimeStampA_cb >= ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB_ax;
            lastU = &_rtDW->LastUAtTimeB_m;
        }
    }
    _rtB->B_97_1931_0 = (_rtB->B_97_1929_0 - *lastU) / (ssGetT(S) - rtb_B_97_47_0);
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    if (_rtB->B_97_1945_0 >= _rtP->P_1861) {
        rtb_B_97_1948_0 = _rtB->B_97_1944_0;
    } else {
        rtb_B_97_1948_0 = rt_Lookup(_rtP->P_1859, 4, ssGetTaskTime(S, 2), _rtP->P_1860);
    }
    if (_rtB->B_97_1949_0 >= _rtP->P_1864) {
        _rtB->B_97_1952_0 = rtb_B_97_1948_0;
    } else {
        _rtB->B_97_1952_0 = _rtB->B_97_1950_0;
    }
}
_rtB->B_97_1954_0 = rt_Lookup(_rtP->P_1865, 5, ssGetT(S), _rtP->P_1866);
if (_rtB->B_97_1943_0 >= _rtP->P_1867) {
    _rtB->B_97_1955_0 = _rtB->B_97_1952_0;
} else {
    _rtB->B_97_1955_0 = _rtB->B_97_1954_0;
}
if ((_rtDW->TimeStampA_co >= ssGetT(S)) && (_rtDW->TimeStampB_f >= ssGetT(S))) {
    _rtB->B_97_1956_0 = 0.0;
} else {
    rtb_B_97_47_0 = _rtDW->TimeStampA_co;
    lastU = &_rtDW->LastUAtTimeA_l;
    if (_rtDW->TimeStampA_co < _rtDW->TimeStampB_f) {
        if (_rtDW->TimeStampB_f < ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB_f;
            lastU = &_rtDW->LastUAtTimeB_i;
        }
    } else {
        if (_rtDW->TimeStampA_co >= ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB_f;
            lastU = &_rtDW->LastUAtTimeB_i;
        }
    }
    _rtB->B_97_1956_0 = (_rtB->B_97_1954_0 - *lastU) / (ssGetT(S) - rtb_B_97_47_0);
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    if (_rtB->B_97_1962_0 >= _rtP->P_1871) {
        _rtB->B_97_1964_0 = rtb_B_97_1948_0;
    } else {
        _rtB->B_97_1964_0 = _rtB->B_97_1950_0;
    }
}
_rtB->B_97_1966_0 = rt_Lookup(_rtP->P_1872, 5, ssGetT(S), _rtP->P_1873);
if (_rtB->B_97_1961_0 >= _rtP->P_1874) {
    _rtB->B_97_1967_0 = _rtB->B_97_1964_0;
} else {
    _rtB->B_97_1967_0 = _rtB->B_97_1966_0;
}
if ((_rtDW->TimeStampA_k >= ssGetT(S)) && (_rtDW->TimeStampB_g >= ssGetT(S))) {
    _rtB->B_97_1968_0 = 0.0;
} else {
    rtb_B_97_47_0 = _rtDW->TimeStampA_k;
    lastU = &_rtDW->LastUAtTimeA_m;
    if (_rtDW->TimeStampA_k < _rtDW->TimeStampB_g) {
        if (_rtDW->TimeStampB_g < ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB_g;
            lastU = &_rtDW->LastUAtTimeB_n;
        }
    } else {
        if (_rtDW->TimeStampA_k >= ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB_g;
            lastU = &_rtDW->LastUAtTimeB_n;
        }
    }
    _rtB->B_97_1968_0 = (_rtB->B_97_1966_0 - *lastU) / (ssGetT(S) - rtb_B_97_47_0);
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    if (_rtB->B_97_1974_0 >= _rtP->P_1878) {
        _rtB->B_97_1976_0 = rtb_B_97_1948_0;
    } else {
        _rtB->B_97_1976_0 = _rtB->B_97_1950_0;
    }
}
_rtB->B_97_1978_0 = rt_Lookup(_rtP->P_1879, 5, ssGetT(S), _rtP->P_1880);
if (_rtB->B_97_1973_0 >= _rtP->P_1881) {
    _rtB->B_97_1979_0 = _rtB->B_97_1976_0;
} else {
    _rtB->B_97_1979_0 = _rtB->B_97_1978_0;
}
if ((_rtDW->TimeStampA_n >= ssGetT(S)) && (_rtDW->TimeStampB_e >= ssGetT(S))) {
    _rtB->B_97_1980_0 = 0.0;
} else {
    rtb_B_97_47_0 = _rtDW->TimeStampA_n;
    lastU = &_rtDW->LastUAtTimeA_n;
    if (_rtDW->TimeStampA_n < _rtDW->TimeStampB_e) {
        if (_rtDW->TimeStampB_e < ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB_e;
            lastU = &_rtDW->LastUAtTimeB_nm;
        }
    } else {
        if (_rtDW->TimeStampA_n >= ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB_e;
            lastU = &_rtDW->LastUAtTimeB_nm;
        }
    }
    _rtB->B_97_1980_0 = (_rtB->B_97_1978_0 - *lastU) / (ssGetT(S) - rtb_B_97_47_0);
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    if (_rtB->B_97_1994_0 >= _rtP->P_1888) {
        rtb_B_97_1997_0 = _rtB->B_97_1993_0;
    } else {
        rtb_B_97_1997_0 = rt_Lookup(_rtP->P_1886, 26, ssGetTaskTime(S, 2), _rtP->P_1887);
    }
    if (_rtB->B_97_1998_0 >= _rtP->P_1891) {
        _rtB->B_97_2001_0 = rtb_B_97_1997_0;
    } else {
        _rtB->B_97_2001_0 = _rtB->B_97_1999_0;
    }
}
_rtB->B_97_2003_0 = rt_Lookup(_rtP->P_1892, 5, ssGetT(S), _rtP->P_1893);
if (_rtB->B_97_1992_0 >= _rtP->P_1894) {
    _rtB->B_97_2004_0 = _rtB->B_97_2001_0;
} else {
    _rtB->B_97_2004_0 = _rtB->B_97_2003_0;
}
if ((_rtDW->TimeStampA_k3 >= ssGetT(S)) && (_rtDW->TimeStampB_j >= ssGetT(S))) {
    _rtB->B_97_2005_0 = 0.0;
} else {
    rtb_B_97_47_0 = _rtDW->TimeStampA_k3;
    lastU = &_rtDW->LastUAtTimeA_mm;
    if (_rtDW->TimeStampA_k3 < _rtDW->TimeStampB_j) {
        if (_rtDW->TimeStampB_j < ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB_j;
            lastU = &_rtDW->LastUAtTimeB_l;
        }
    } else {
        if (_rtDW->TimeStampA_k3 >= ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB_j;
            lastU = &_rtDW->LastUAtTimeB_l;
        }
    }
    _rtB->B_97_2005_0 = (_rtB->B_97_2003_0 - *lastU) / (ssGetT(S) - rtb_B_97_47_0);
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    if (_rtB->B_97_2011_0 >= _rtP->P_1898) {
        _rtB->B_97_2013_0 = rtb_B_97_1997_0;
    } else {
        _rtB->B_97_2013_0 = _rtB->B_97_1999_0;
    }
}
_rtB->B_97_2015_0 = rt_Lookup(_rtP->P_1899, 5, ssGetT(S), _rtP->P_1900);
if (_rtB->B_97_2010_0 >= _rtP->P_1901) {
    _rtB->B_97_2016_0 = _rtB->B_97_2013_0;
} else {
    _rtB->B_97_2016_0 = _rtB->B_97_2015_0;
}
if ((_rtDW->TimeStampA_l >= ssGetT(S)) && (_rtDW->TimeStampB_l >= ssGetT(S))) {
    _rtB->B_97_2017_0 = 0.0;
} else {
    rtb_B_97_47_0 = _rtDW->TimeStampA_l;
    lastU = &_rtDW->LastUAtTimeA_h;
    if (_rtDW->TimeStampA_l < _rtDW->TimeStampB_l) {
        if (_rtDW->TimeStampB_l < ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB_l;
            lastU = &_rtDW->LastUAtTimeB_b;
        }
    } else {
        if (_rtDW->TimeStampA_l >= ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB_l;
            lastU = &_rtDW->LastUAtTimeB_b;
        }
    }
    _rtB->B_97_2017_0 = (_rtB->B_97_2015_0 - *lastU) / (ssGetT(S) - rtb_B_97_47_0);
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    if (_rtB->B_97_2023_0 >= _rtP->P_1905) {
        _rtB->B_97_2025_0 = rtb_B_97_1997_0;
    } else {
        _rtB->B_97_2025_0 = _rtB->B_97_1999_0;
    }
}
_rtB->B_97_2027_0 = rt_Lookup(_rtP->P_1906, 5, ssGetT(S), _rtP->P_1907);
if (_rtB->B_97_2022_0 >= _rtP->P_1908) {
    _rtB->B_97_2028_0 = _rtB->B_97_2025_0;
} else {
    _rtB->B_97_2028_0 = _rtB->B_97_2027_0;
}
if ((_rtDW->TimeStampA_cm >= ssGetT(S)) && (_rtDW->TimeStampB_gf >= ssGetT(S))) {
    _rtB->B_97_2029_0 = 0.0;
} else {
    rtb_B_97_47_0 = _rtDW->TimeStampA_cm;
    lastU = &_rtDW->LastUAtTimeA_n3;
    if (_rtDW->TimeStampA_cm < _rtDW->TimeStampB_gf) {
        if (_rtDW->TimeStampB_gf < ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB_gf;
            lastU = &_rtDW->LastUAtTimeB_g;
        }
    } else {
        if (_rtDW->TimeStampA_cm >= ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB_gf;
            lastU = &_rtDW->LastUAtTimeB_g;
        }
    }
    _rtB->B_97_2029_0 = (_rtB->B_97_2027_0 - *lastU) / (ssGetT(S) - rtb_B_97_47_0);
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    if (_rtB->B_97_2043_0 >= _rtP->P_1915) {
        rtb_B_97_2046_0 = _rtB->B_97_2042_0;
    } else {
        rtb_B_97_2046_0 = rt_Lookup(_rtP->P_1913, 4, ssGetTaskTime(S, 2), _rtP->P_1914);
    }
    if (_rtB->B_97_2047_0 >= _rtP->P_1918) {
        _rtB->B_97_2050_0 = rtb_B_97_2046_0;
    } else {
        _rtB->B_97_2050_0 = _rtB->B_97_2048_0;
    }
}
_rtB->B_97_2052_0 = rt_Lookup(_rtP->P_1919, 5, ssGetT(S), _rtP->P_1920);
if (_rtB->B_97_2041_0 >= _rtP->P_1921) {
    _rtB->B_97_2053_0 = _rtB->B_97_2050_0;
} else {
    _rtB->B_97_2053_0 = _rtB->B_97_2052_0;
}
if ((_rtDW->TimeStampA_f >= ssGetT(S)) && (_rtDW->TimeStampB_gp >= ssGetT(S))) {
    _rtB->B_97_2054_0 = 0.0;
} else {
    rtb_B_97_47_0 = _rtDW->TimeStampA_f;
    lastU = &_rtDW->LastUAtTimeA_dj;
    if (_rtDW->TimeStampA_f < _rtDW->TimeStampB_gp) {
        if (_rtDW->TimeStampB_gp < ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB_gp;
            lastU = &_rtDW->LastUAtTimeB_o;
        }
    } else {
        if (_rtDW->TimeStampA_f >= ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB_gp;
            lastU = &_rtDW->LastUAtTimeB_o;
        }
    }
    _rtB->B_97_2054_0 = (_rtB->B_97_2052_0 - *lastU) / (ssGetT(S) - rtb_B_97_47_0);
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    if (_rtB->B_97_2060_0 >= _rtP->P_1925) {
        _rtB->B_97_2062_0 = rtb_B_97_2046_0;
    } else {
        _rtB->B_97_2062_0 = _rtB->B_97_2048_0;
    }
}
_rtB->B_97_2064_0 = rt_Lookup(_rtP->P_1926, 5, ssGetT(S), _rtP->P_1927);
if (_rtB->B_97_2059_0 >= _rtP->P_1928) {
    _rtB->B_97_2065_0 = _rtB->B_97_2062_0;
} else {
    _rtB->B_97_2065_0 = _rtB->B_97_2064_0;
}
if ((_rtDW->TimeStampA_d >= ssGetT(S)) && (_rtDW->TimeStampB_eg >= ssGetT(S))) {
    _rtB->B_97_2066_0 = 0.0;
} else {
    rtb_B_97_47_0 = _rtDW->TimeStampA_d;
    lastU = &_rtDW->LastUAtTimeA_dy;
    if (_rtDW->TimeStampA_d < _rtDW->TimeStampB_eg) {
        if (_rtDW->TimeStampB_eg < ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB_eg;
            lastU = &_rtDW->LastUAtTimeB_a;
        }
    } else {
        if (_rtDW->TimeStampA_d >= ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB_eg;
            lastU = &_rtDW->LastUAtTimeB_a;
        }
    }
    _rtB->B_97_2066_0 = (_rtB->B_97_2064_0 - *lastU) / (ssGetT(S) - rtb_B_97_47_0);
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    if (_rtB->B_97_2072_0 >= _rtP->P_1932) {
        _rtB->B_97_2074_0 = rtb_B_97_2046_0;
    } else {
        _rtB->B_97_2074_0 = _rtB->B_97_2048_0;
    }
}
_rtB->B_97_2076_0 = rt_Lookup(_rtP->P_1933, 5, ssGetT(S), _rtP->P_1934);
if (_rtB->B_97_2071_0 >= _rtP->P_1935) {
    _rtB->B_97_2077_0 = _rtB->B_97_2074_0;
} else {
    _rtB->B_97_2077_0 = _rtB->B_97_2076_0;
}
if ((_rtDW->TimeStampA_b >= ssGetT(S)) && (_rtDW->TimeStampB_k >= ssGetT(S))) {
    _rtB->B_97_2078_0 = 0.0;
} else {
    rtb_B_97_47_0 = _rtDW->TimeStampA_b;
    lastU = &_rtDW->LastUAtTimeA_f;
    if (_rtDW->TimeStampA_b < _rtDW->TimeStampB_k) {
        if (_rtDW->TimeStampB_k < ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB_k;
            lastU = &_rtDW->LastUAtTimeB_h;
        }
    } else {
        if (_rtDW->TimeStampA_b >= ssGetT(S)) {
            rtb_B_97_47_0 = _rtDW->TimeStampB_k;
            lastU = &_rtDW->LastUAtTimeB_h;
        }
    }
    _rtB->B_97_2078_0 = (_rtB->B_97_2076_0 - *lastU) / (ssGetT(S) - rtb_B_97_47_0);
}








    
  

                    UNUSED_PARAMETER(tid);



  
  
    
        }
      
  






    
  

          /* Outputs for root system: '<Root>' */
      
            
  
  
      
  
     static  void mdlOutputsTID5(SimStruct *S, int_T tid)
  {
  
              /* local block i/o variables */
    
                real_T B_97_155_0;



  
    
int32_T i;
B_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *_rtB;
P_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *_rtP;
DW_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *_rtDW;

    
  

          
  
    
  
    
  


                   
                  

_rtDW = ((DW_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetRootDWork(S));
_rtP = ((P_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetModelRtp(S));
_rtB = ((B_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) _ssGetModelBlockIO(S));
memcpy(&_rtB->B_97_0_0[0], &_rtP->P_537[0], 13U * sizeof(real_T));
_rtB->B_97_7_0 = _rtP->P_542;


      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_SaturationTID5(S, &_rtB->Saturation_i, &_rtP->Saturation_i);
_rtB->B_97_10_0 = _rtP->P_543;
_rtB->B_97_15_0 = _rtP->P_548;
_rtB->B_97_18_0 = _rtP->P_551;
for (i = 0; i < 25; i++) {
    _rtB->B_97_11_0[i] = _rtP->P_544[i];
    _rtB->B_97_13_0[i] = _rtP->P_546[i];
    _rtB->B_97_14_0[i] = _rtP->P_547[i];
    _rtB->B_97_19_0[i] = _rtP->P_552[i];
    _rtB->B_97_20_0[i] = _rtP->P_553 * _rtB->B_97_13_0[i];
}
_rtB->B_97_38_0 = _rtP->P_562;


      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_SaturationTID5(S, &_rtB->Saturation, &_rtP->Saturation);
_rtB->B_97_41_0 = _rtP->P_563;
_rtB->B_97_46_0 = _rtP->P_568;
_rtB->B_97_49_0 = _rtP->P_571;
for (i = 0; i < 25; i++) {
    _rtB->B_97_42_0[i] = _rtP->P_564[i];
    _rtB->B_97_44_0[i] = _rtP->P_566[i];
    _rtB->B_97_45_0[i] = _rtP->P_567[i];
    _rtB->B_97_50_0[i] = _rtP->P_572[i];
    _rtB->B_97_51_0[i] = _rtP->P_573 * _rtB->B_97_44_0[i];
}
_rtB->B_97_63_0 = _rtP->P_578;
_rtB->B_97_64_0 = _rtP->P_579;
_rtB->B_97_65_0 = _rtP->P_580;
_rtB->B_97_77_0 = _rtP->P_592;
_rtB->B_97_85_0 = _rtP->P_594;
_rtB->B_97_86_0 = _rtP->P_595;
_rtB->B_97_96_0 = _rtP->P_598;
_rtB->B_97_97_0 = _rtP->P_599;
_rtB->B_97_98_0 = _rtP->P_600;
_rtB->B_97_104_0 = _rtP->P_601;
_rtB->B_97_135_0 = _rtP->P_619;
B_97_155_0 = _rtP->P_628;
_rtB->B_97_159_0 = _rtP->P_629;
_rtB->B_97_161_0 = _rtP->P_630;
_rtB->B_97_178_0 = _rtP->P_637;
_rtB->B_97_180_0 = _rtP->P_638;
_rtB->B_97_217_0 = _rtP->P_664;
_rtB->B_97_219_0 = _rtP->P_665;
_rtB->B_97_352_0 = _rtP->P_805;
_rtB->B_97_355_0 = _rtP->P_806;
_rtB->B_97_358_0 = _rtP->P_807;
for (i = 0; i < 6; i++) {
    _rtB->B_97_372_0[i] = _rtP->P_817[i];
    _rtB->B_97_373_0[i] = _rtP->P_818[i];
}
_rtB->B_97_379_0 = _rtP->P_820;
_rtB->B_97_398_0 = _rtP->P_825;
_rtB->B_97_399_0 = _rtP->P_826;
_rtB->B_97_400_0 = _rtP->P_827;
_rtB->B_97_402_0 = _rtP->P_828;
_rtB->B_27_1_0 = _rtP->P_60;
_rtB->B_27_4_0 = _rtP->P_61;
_rtB->B_27_34_0 = _rtP->P_84;
_rtB->B_27_36_0 = _rtP->P_86;
_rtB->B_27_46_0 = _rtP->P_95;
_rtB->B_27_49_0 = _rtP->P_96;
_rtB->B_27_77_0 = _rtP->P_119;
_rtB->B_27_82_0 = _rtP->P_121;
_rtB->B_27_84_0 = _rtP->P_122;
_rtB->B_27_96_0 = _rtP->P_125;
_rtB->B_27_98_0 = _rtP->P_126;
_rtB->B_27_140_0 = _rtP->P_154 * _rtP->P_153 / _rtP->P_155;
_rtB->B_27_145_0[0] = _rtP->P_158[0];
_rtB->B_27_145_0[1] = _rtP->P_158[1];
_rtB->B_27_145_0[2] = _rtP->P_158[2];
_rtB->B_27_146_0 = _rtP->P_159;
_rtB->B_27_147_0 = _rtP->P_160;
_rtB->B_27_165_0 = _rtP->P_162;
_rtB->B_27_167_0 = _rtP->P_163;
_rtB->B_27_188_0 = _rtP->P_176;
_rtB->B_23_0_0 = _rtP->P_27;
_rtB->B_23_2_0 = _rtP->P_28;
if (ssIsMajorTimeStep(S) != 0) {
    srUpdateBC(_rtDW->AutomaticGainControl_SubsysRanBC_nm);
}
_rtB->B_27_240_0 = _rtP->P_216;
_rtB->B_27_242_0 = _rtP->P_218;
if (ssIsMajorTimeStep(S) != 0) {
    srUpdateBC(_rtDW->GridSupportingasCurrentSourceGridFeeding_SubsysRanBC);
}
_rtB->B_97_413_0 = _rtP->P_836;
_rtB->B_97_414_0 = _rtP->P_837;
_rtB->B_36_0_0 = _rtP->P_256;
_rtB->B_36_9_0 = _rtP->P_264;
_rtB->B_36_12_0 = _rtP->P_266 * _rtP->P_265;
_rtB->B_36_16_0 = _rtP->P_269;
_rtB->B_36_18_0 = _rtP->P_270;
_rtB->B_36_25_0 = _rtP->P_278;
_rtB->B_36_31_0 = (uint8_T)(_rtP->P_283 == _rtP->P_284);
_rtB->B_36_35_0 = (uint8_T)(_rtP->P_283 == _rtP->P_285);
_rtB->B_36_51_0 = _rtP->P_293;
_rtB->B_36_55_0 = _rtP->P_294;
_rtB->B_36_57_0 = _rtP->P_295;
_rtB->B_36_78_0 = _rtP->P_308;
_rtB->B_30_0_0 = _rtP->P_223;
_rtB->B_30_2_0 = _rtP->P_224;
if (ssIsMajorTimeStep(S) != 0) {
    srUpdateBC(_rtDW->AutomaticGainControl_SubsysRanBC_p);
}
_rtB->B_36_107_0 = (uint8_T)(_rtP->P_333 == _rtP->P_334);
_rtB->B_36_111_0 = (uint8_T)(_rtP->P_333 == _rtP->P_335);
if (ssIsMajorTimeStep(S) != 0) {
    srUpdateBC(_rtDW->Gridforming_SubsysRanBC);
}
_rtB->B_97_424_0 = _rtP->P_838;
_rtB->B_97_442_0 = _rtP->P_856;
_rtB->B_97_446_0 = _rtP->P_857;
_rtB->B_97_448_0 = _rtP->P_858;
_rtB->B_97_469_0 = _rtP->P_871;
_rtB->B_44_0_0 = _rtP->P_340;
_rtB->B_44_2_0 = _rtP->P_341;
if (ssIsMajorTimeStep(S) != 0) {
    srUpdateBC(_rtDW->AutomaticGainControl_SubsysRanBC_hn);
}
_rtB->B_97_576_0 = _rtP->P_908;
_rtB->B_97_578_0 = _rtP->P_909;
_rtB->B_97_615_0 = _rtP->P_935;
_rtB->B_97_617_0 = _rtP->P_936;
_rtB->B_97_750_0 = _rtP->P_1076;
_rtB->B_97_753_0 = _rtP->P_1077;
_rtB->B_97_756_0 = _rtP->P_1078;
for (i = 0; i < 6; i++) {
    _rtB->B_97_770_0[i] = _rtP->P_1088[i];
    _rtB->B_97_771_0[i] = _rtP->P_1089[i];
}
_rtB->B_97_788_0 = _rtP->P_1098;
_rtB->B_97_810_0 = _rtP->P_1103;
_rtB->B_97_828_0 = _rtP->P_1121;
_rtB->B_97_832_0 = _rtP->P_1122;
_rtB->B_97_834_0 = _rtP->P_1123;
_rtB->B_97_855_0 = _rtP->P_1136;
_rtB->B_66_0_0 = _rtP->P_420;
_rtB->B_66_2_0 = _rtP->P_421;
if (ssIsMajorTimeStep(S) != 0) {
    srUpdateBC(_rtDW->AutomaticGainControl_SubsysRanBC_n);
}
_rtB->B_97_890_0 = _rtP->P_1167;
_rtB->B_97_893_0 = _rtP->P_1168;
_rtB->B_97_921_0 = _rtP->P_1190;
_rtB->B_97_924_0 = _rtP->P_1191;
_rtB->B_97_952_0 = _rtP->P_1210;
_rtB->B_97_954_0 = _rtP->P_1212;
_rtB->B_97_956_0 = _rtP->P_1214;
_rtB->B_97_959_0 = _rtP->P_1217;
_rtB->B_97_963_0 = _rtP->P_1222;
_rtB->B_97_968_0 = _rtP->P_1224;
_rtB->B_97_970_0 = _rtP->P_1225;
_rtB->B_97_982_0 = _rtP->P_1228;
_rtB->B_97_984_0 = _rtP->P_1229;
_rtB->B_97_1026_0 = _rtP->P_1257 * _rtP->P_1256 / _rtP->P_1258;
_rtB->B_97_1031_0[0] = _rtP->P_1261[0];
_rtB->B_97_1031_0[1] = _rtP->P_1261[1];
_rtB->B_97_1031_0[2] = _rtP->P_1261[2];
_rtB->B_97_1032_0 = _rtP->P_1262;
_rtB->B_97_1033_0 = _rtP->P_1263;
_rtB->B_97_1050_0 = _rtP->P_1265;
_rtB->B_97_1052_0 = _rtP->P_1266;
_rtB->B_97_1073_0 = _rtP->P_1279;
_rtB->B_93_0_0 = _rtP->P_506;
_rtB->B_93_2_0 = _rtP->P_507;
if (ssIsMajorTimeStep(S) != 0) {
    srUpdateBC(_rtDW->AutomaticGainControl_SubsysRanBC);
}
_rtB->B_97_1125_0 = _rtP->P_1319;
_rtB->B_97_1127_0 = _rtP->P_1321;
_rtB->B_97_1163_0 = _rtP->P_1323;
_rtB->B_97_1164_0 = _rtP->P_1324;
_rtB->B_97_1165_0 = _rtP->P_1937;
_rtB->B_97_1166_0 = _rtP->P_1325;
_rtB->B_97_1169_0 = _rtP->P_1327;
_rtB->B_97_1181_0 = _rtP->P_1332;
_rtB->B_97_1170_0[0] = _rtP->P_1328[0];
_rtB->B_97_1182_0[0] = _rtP->P_1333[0];
_rtB->B_97_1170_0[1] = _rtP->P_1328[1];
_rtB->B_97_1182_0[1] = _rtP->P_1333[1];
_rtB->B_97_1203_0 = _rtP->P_1356;
_rtB->B_97_1210_0 = (real_T)(_rtP->P_1358 == _rtP->P_1359) * _rtP->P_1360;
_rtB->B_97_1245_0 = _rtP->P_1389;
_rtB->B_97_1251_0[0] = _rtP->P_1396[0];
_rtB->B_97_1251_0[1] = _rtP->P_1396[1];
_rtB->B_97_1355_0 = _rtP->P_1516;
_rtB->B_97_1358_0 = _rtP->P_1517;
_rtB->B_97_1361_0 = _rtP->P_1518;
_rtB->B_97_1373_0 = _rtP->P_1522;
_rtB->B_97_1391_0 = _rtP->P_1540;
_rtB->B_97_1395_0 = _rtP->P_1541;
_rtB->B_97_1397_0 = _rtP->P_1542;
_rtB->B_97_1418_0 = _rtP->P_1555;
_rtB->B_55_0_0 = _rtP->P_380;
_rtB->B_55_2_0 = _rtP->P_381;
if (ssIsMajorTimeStep(S) != 0) {
    srUpdateBC(_rtDW->AutomaticGainControl_SubsysRanBC_e);
}
_rtB->B_97_1566_0 = _rtP->P_1691;
_rtB->B_97_1569_0 = _rtP->P_1692;
_rtB->B_97_1572_0 = _rtP->P_1693;
_rtB->B_97_1584_0 = _rtP->P_1696;
_rtB->B_97_1602_0 = _rtP->P_1714;
_rtB->B_97_1606_0 = _rtP->P_1715;
_rtB->B_97_1608_0 = _rtP->P_1716;
_rtB->B_97_1628_0 = _rtP->P_1728;
_rtB->B_77_0_0 = _rtP->P_460;
_rtB->B_77_2_0 = _rtP->P_461;
if (ssIsMajorTimeStep(S) != 0) {
    srUpdateBC(_rtDW->AutomaticGainControl_SubsysRanBC_h);
}
_rtB->B_97_1728_0 = _rtP->P_1938;
_rtB->B_97_1729_0 = _rtP->P_1760;
_rtB->B_97_1732_0 = _rtP->P_1762;
_rtB->B_97_1744_0 = _rtP->P_1767;
_rtB->B_97_1733_0[0] = _rtP->P_1763[0];
_rtB->B_97_1745_0[0] = _rtP->P_1768[0];
_rtB->B_97_1733_0[1] = _rtP->P_1763[1];
_rtB->B_97_1745_0[1] = _rtP->P_1768[1];
_rtB->B_97_1761_0 = _rtP->P_1774;
_rtB->B_97_1762_0 = _rtP->P_1775;
_rtB->B_97_1770_0 = _rtP->P_1793;
_rtB->B_97_1777_0 = (real_T)(_rtP->P_1795 == _rtP->P_1796) * _rtP->P_1797;
_rtB->B_97_1797_0 = _rtP->P_1816;
_rtB->B_97_1803_0[0] = _rtP->P_1823[0];
_rtB->B_97_1803_0[1] = _rtP->P_1823[1];
_rtB->B_97_1894_0 = _rtP->P_1829;
_rtB->B_97_1895_0 = _rtP->P_1830;
_rtB->B_97_1896_0 = _rtP->P_1831;
_rtB->B_97_1900_0 = _rtP->P_1835;
_rtB->B_97_1901_0 = _rtP->P_1836;
_rtB->B_97_1912_0 = _rtP->P_1842;
_rtB->B_97_1913_0 = _rtP->P_1843;
_rtB->B_97_1924_0 = _rtP->P_1849;
_rtB->B_97_1925_0 = _rtP->P_1850;
_rtB->B_97_1943_0 = _rtP->P_1856;
_rtB->B_97_1944_0 = _rtP->P_1857;
_rtB->B_97_1945_0 = _rtP->P_1858;
_rtB->B_97_1949_0 = _rtP->P_1862;
_rtB->B_97_1950_0 = _rtP->P_1863;
_rtB->B_97_1961_0 = _rtP->P_1869;
_rtB->B_97_1962_0 = _rtP->P_1870;
_rtB->B_97_1973_0 = _rtP->P_1876;
_rtB->B_97_1974_0 = _rtP->P_1877;
_rtB->B_97_1992_0 = _rtP->P_1883;
_rtB->B_97_1993_0 = _rtP->P_1884;
_rtB->B_97_1994_0 = _rtP->P_1885;
_rtB->B_97_1998_0 = _rtP->P_1889;
_rtB->B_97_1999_0 = _rtP->P_1890;
_rtB->B_97_2010_0 = _rtP->P_1896;
_rtB->B_97_2011_0 = _rtP->P_1897;
_rtB->B_97_2022_0 = _rtP->P_1903;
_rtB->B_97_2023_0 = _rtP->P_1904;
_rtB->B_97_2041_0 = _rtP->P_1910;
_rtB->B_97_2042_0 = _rtP->P_1911;
_rtB->B_97_2043_0 = _rtP->P_1912;
_rtB->B_97_2047_0 = _rtP->P_1916;
_rtB->B_97_2048_0 = _rtP->P_1917;
_rtB->B_97_2059_0 = _rtP->P_1923;
_rtB->B_97_2060_0 = _rtP->P_1924;
_rtB->B_97_2071_0 = _rtP->P_1930;
_rtB->B_97_2072_0 = _rtP->P_1931;






    
  

                UNUSED_PARAMETER(tid);



  
  
    
        }
      
  








       
    
  

          /* Update for root system: '<Root>' */
      
            #define MDL_UPDATE
    
  
  
      
  
     static  void mdlUpdate(SimStruct *S, int_T tid)
  {
  
      
        
    
real_T HoldSine;
real_T *lastU;
int32_T i;
B_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *_rtB;
P_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *_rtP;
X_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *_rtX;
DW_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *_rtDW;
  
    
  

                
      
    
  

                      
                  
_rtDW = ((DW_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetRootDWork(S));
_rtX = ((X_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetContStates(S));
_rtP = ((P_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetModelRtp(S));
_rtB = ((B_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) _ssGetModelBlockIO(S));
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtDW->DiscreteTimeIntegrator1_DSTATE += _rtP->P_538 * _rtB->B_97_1258_0;
    for (i = 0; i < 5; i++) {
        _rtDW->fluxes_DSTATE[i] = _rtB->B_97_24_0[i];
    }


      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_Saturation_Update(S, _rtB->B_97_8_0, &_rtB->Saturation_i, &_rtDW->Saturation_i);
    _rtDW->DiscreteTimeIntegrator_DSTATE += _rtP->P_549 * _rtB->B_97_1257_0;
    _rtDW->voltages_DSTATE[0] = _rtB->B_97_1179_0;
    _rtDW->voltages_DSTATE[1] = _rtB->B_97_1180_0;
    _rtDW->voltages_DSTATE[2] = _rtB->B_97_1250_0;
    _rtDW->voltages_DSTATE[3] = _rtB->B_97_1251_0[0];
    _rtDW->voltages_DSTATE[4] = _rtB->B_97_1251_0[1];
    _rtDW->DiscreteTimeIntegrator1_DSTATE_h += _rtP->P_558 * _rtB->B_97_1810_0;
    for (i = 0; i < 5; i++) {
        _rtDW->fluxes_DSTATE_i[i] = _rtB->B_97_55_0[i];
    }


      Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_Saturation_Update(S, _rtB->B_97_39_0, &_rtB->Saturation, &_rtDW->Saturation);
    _rtDW->DiscreteTimeIntegrator_DSTATE_f += _rtP->P_569 * _rtB->B_97_1809_0;
    _rtDW->voltages_DSTATE_c[0] = _rtB->B_97_1742_0;
    _rtDW->voltages_DSTATE_c[1] = _rtB->B_97_1743_0;
    _rtDW->voltages_DSTATE_c[2] = _rtB->B_97_1802_0;
    _rtDW->voltages_DSTATE_c[3] = _rtB->B_97_1803_0[0];
    _rtDW->voltages_DSTATE_c[4] = _rtB->B_97_1803_0[1];
}
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    _rtDW->itinit1_PreviousInput = _rtB->B_97_135_0;
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtDW->Currentfilter_states = (_rtB->B_97_129_0 - _rtP->P_584[1] * _rtDW->Currentfilter_states) / _rtP->P_584[0];
}
i = ssIsSampleHit(S, 1, 0);
if (i != 0) {
    _rtDW->itinit_PreviousInput = _rtB->B_97_146_0;
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtDW->inti_IC_LOADING = 0U;
    _rtDW->inti_DSTATE += _rtP->P_588 * _rtB->B_97_129_0;
    if (_rtDW->inti_DSTATE >= _rtP->P_589) {
        _rtDW->inti_DSTATE = _rtP->P_589;
    } else {
        if (_rtDW->inti_DSTATE <= _rtP->P_590) {
            _rtDW->inti_DSTATE = _rtP->P_590;
        }
    }
    if (_rtB->B_97_72_0 > 0.0) {
        _rtDW->inti_PrevResetState = 1;
    } else if (_rtB->B_97_72_0 < 0.0) {
        _rtDW->inti_PrevResetState = -1;
    } else if (_rtB->B_97_72_0 == 0.0) {
        _rtDW->inti_PrevResetState = 0;
    } else {
        _rtDW->inti_PrevResetState = 2;
    }
    _rtDW->DiscreteTimeIntegrator_DSTATE_l += _rtP->P_602 * _rtB->B_97_143_0;
    _rtDW->Memory2_PreviousInput = _rtB->B_97_147_0;
                    /* Level2 S-Function Block: '<S1410>/B_82_9' (sfun_spssw_discc) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 116, SS_CALL_MDL_UPDATE);

     _rtDW->Integ4_SYSTEM_ENABLE = 0U;
    _rtDW->Integ4_DSTATE = _rtP->P_639 * _rtB->B_97_188_0 + _rtB->B_97_189_0;
                    /* Level2 S-Function Block: '<S168>/B_82_31' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 191, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE = _rtB->B_97_188_0;
    _rtDW->Integ4_SYSTEM_ENABLE_a = 0U;
    _rtDW->Integ4_DSTATE_f = _rtP->P_649 * _rtB->B_97_198_0 + _rtB->B_97_199_0;
                    /* Level2 S-Function Block: '<S169>/B_82_33' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 201, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_d = _rtB->B_97_198_0;
    _rtDW->Integ4_SYSTEM_ENABLE_c = 0U;
    _rtDW->Integ4_DSTATE_m = _rtP->P_666 * _rtB->B_97_227_0 + _rtB->B_97_228_0;
                    /* Level2 S-Function Block: '<S162>/B_82_35' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 230, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_e = _rtB->B_97_227_0;
    _rtDW->Integ4_SYSTEM_ENABLE_l = 0U;
    _rtDW->Integ4_DSTATE_d = _rtP->P_676 * _rtB->B_97_237_0 + _rtB->B_97_238_0;
                    /* Level2 S-Function Block: '<S163>/B_82_37' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 240, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_f = _rtB->B_97_237_0;
    HoldSine = _rtDW->lastSin;
    _rtDW->lastSin = _rtDW->lastSin * _rtP->P_699 + _rtDW->lastCos * _rtP->P_698;
    _rtDW->lastCos = _rtDW->lastCos * _rtP->P_699 - HoldSine * _rtP->P_698;
    _rtDW->Integ4_SYSTEM_ENABLE_p = 0U;
    _rtDW->Integ4_DSTATE_k[0] = _rtP->P_702 * _rtB->B_97_267_0[0] + _rtB->B_97_268_0[0];
    _rtDW->Integ4_DSTATE_k[1] = _rtP->P_702 * _rtB->B_97_267_0[1] + _rtB->B_97_268_0[1];
    _rtDW->Integ4_DSTATE_k[2] = _rtP->P_702 * _rtB->B_97_267_0[2] + _rtB->B_97_268_0[2];
                    /* Level2 S-Function Block: '<S694>/B_82_39' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 270, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_dz[0] = _rtB->B_97_267_0[0];
    _rtDW->UnitDelay1_DSTATE[0] = _rtB->B_97_277_0[0];
    _rtDW->UnitDelay_DSTATE_dz[1] = _rtB->B_97_267_0[1];
    _rtDW->UnitDelay1_DSTATE[1] = _rtB->B_97_277_0[1];
    _rtDW->UnitDelay_DSTATE_dz[2] = _rtB->B_97_267_0[2];
    _rtDW->UnitDelay1_DSTATE[2] = _rtB->B_97_277_0[2];
    HoldSine = _rtDW->lastSin_m;
    _rtDW->lastSin_m = _rtDW->lastSin_m * _rtP->P_720 + _rtDW->lastCos_i * _rtP->P_719;
    _rtDW->lastCos_i = _rtDW->lastCos_i * _rtP->P_720 - HoldSine * _rtP->P_719;
    _rtDW->Integ4_SYSTEM_ENABLE_f = 0U;
    _rtDW->Integ4_DSTATE_h[0] = _rtP->P_723 * _rtB->B_97_279_0[0] + _rtB->B_97_280_0[0];
    _rtDW->Integ4_DSTATE_h[1] = _rtP->P_723 * _rtB->B_97_279_0[1] + _rtB->B_97_280_0[1];
    _rtDW->Integ4_DSTATE_h[2] = _rtP->P_723 * _rtB->B_97_279_0[2] + _rtB->B_97_280_0[2];
                    /* Level2 S-Function Block: '<S692>/B_82_41' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 282, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_o[0] = _rtB->B_97_279_0[0];
    _rtDW->UnitDelay1_DSTATE_c[0] = _rtB->B_97_289_0[0];
    _rtDW->UnitDelay_DSTATE_o[1] = _rtB->B_97_279_0[1];
    _rtDW->UnitDelay1_DSTATE_c[1] = _rtB->B_97_289_0[1];
    _rtDW->UnitDelay_DSTATE_o[2] = _rtB->B_97_279_0[2];
    _rtDW->UnitDelay1_DSTATE_c[2] = _rtB->B_97_289_0[2];
    HoldSine = _rtDW->lastSin_e;
    _rtDW->lastSin_e = _rtDW->lastSin_e * _rtP->P_745 + _rtDW->lastCos_j * _rtP->P_744;
    _rtDW->lastCos_j = _rtDW->lastCos_j * _rtP->P_745 - HoldSine * _rtP->P_744;
    _rtDW->Integ4_SYSTEM_ENABLE_h = 0U;
    _rtDW->Integ4_DSTATE_n[0] = _rtP->P_748 * _rtB->B_97_297_0[0] + _rtB->B_97_298_0[0];
    _rtDW->Integ4_DSTATE_n[1] = _rtP->P_748 * _rtB->B_97_297_0[1] + _rtB->B_97_298_0[1];
    _rtDW->Integ4_DSTATE_n[2] = _rtP->P_748 * _rtB->B_97_297_0[2] + _rtB->B_97_298_0[2];
                    /* Level2 S-Function Block: '<S700>/B_82_43' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 300, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_ee[0] = _rtB->B_97_297_0[0];
    _rtDW->UnitDelay1_DSTATE_h[0] = _rtB->B_97_307_0[0];
    _rtDW->UnitDelay_DSTATE_ee[1] = _rtB->B_97_297_0[1];
    _rtDW->UnitDelay1_DSTATE_h[1] = _rtB->B_97_307_0[1];
    _rtDW->UnitDelay_DSTATE_ee[2] = _rtB->B_97_297_0[2];
    _rtDW->UnitDelay1_DSTATE_h[2] = _rtB->B_97_307_0[2];
    HoldSine = _rtDW->lastSin_mo;
    _rtDW->lastSin_mo = _rtDW->lastSin_mo * _rtP->P_766 + _rtDW->lastCos_f * _rtP->P_765;
    _rtDW->lastCos_f = _rtDW->lastCos_f * _rtP->P_766 - HoldSine * _rtP->P_765;
    _rtDW->Integ4_SYSTEM_ENABLE_k = 0U;
    _rtDW->Integ4_DSTATE_l[0] = _rtP->P_769 * _rtB->B_97_309_0[0] + _rtB->B_97_310_0[0];
    _rtDW->Integ4_DSTATE_l[1] = _rtP->P_769 * _rtB->B_97_309_0[1] + _rtB->B_97_310_0[1];
    _rtDW->Integ4_DSTATE_l[2] = _rtP->P_769 * _rtB->B_97_309_0[2] + _rtB->B_97_310_0[2];
                    /* Level2 S-Function Block: '<S698>/B_82_45' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 312, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_p[0] = _rtB->B_97_309_0[0];
    _rtDW->UnitDelay1_DSTATE_j[0] = _rtB->B_97_319_0[0];
    _rtDW->UnitDelay_DSTATE_p[1] = _rtB->B_97_309_0[1];
    _rtDW->UnitDelay1_DSTATE_j[1] = _rtB->B_97_319_0[1];
    _rtDW->UnitDelay_DSTATE_p[2] = _rtB->B_97_309_0[2];
    _rtDW->UnitDelay1_DSTATE_j[2] = _rtB->B_97_319_0[2];
}
            {
    real_T **uBuffer = (real_T**)&_rtDW->T_PWORK.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T_PWORK.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T_IWORK.Head = ((_rtDW->T_IWORK.Head < (_rtDW->T_IWORK.CircularBufSize-1)) ? (_rtDW->T_IWORK.Head+1) : 0);

      if (_rtDW->T_IWORK.Head == _rtDW->T_IWORK.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T_IWORK.CircularBufSize, &_rtDW->T_IWORK.Tail, &_rtDW->T_IWORK.Head, &_rtDW->T_IWORK.Last, simTime - _rtP->P_788, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T_IWORK.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T_IWORK.Head] = simTime;
        (*uBuffer)[_rtDW->T_IWORK.Head] = _rtB->B_97_331_0;    
  }
             {
    real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T1_IWORK.Head = ((_rtDW->T1_IWORK.Head < (_rtDW->T1_IWORK.CircularBufSize-1)) ? (_rtDW->T1_IWORK.Head+1) : 0);

      if (_rtDW->T1_IWORK.Head == _rtDW->T1_IWORK.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T1_IWORK.CircularBufSize, &_rtDW->T1_IWORK.Tail, &_rtDW->T1_IWORK.Head, &_rtDW->T1_IWORK.Last, simTime - _rtP->P_791, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T1_IWORK.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T1_IWORK.Head] = simTime;
        (*uBuffer)[_rtDW->T1_IWORK.Head] = _rtB->B_97_334_0;    
  }
             {
    real_T **uBuffer = (real_T**)&_rtDW->T_PWORK_l.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T_PWORK_l.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T_IWORK_h.Head = ((_rtDW->T_IWORK_h.Head < (_rtDW->T_IWORK_h.CircularBufSize-1)) ? (_rtDW->T_IWORK_h.Head+1) : 0);

      if (_rtDW->T_IWORK_h.Head == _rtDW->T_IWORK_h.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T_IWORK_h.CircularBufSize, &_rtDW->T_IWORK_h.Tail, &_rtDW->T_IWORK_h.Head, &_rtDW->T_IWORK_h.Last, simTime - _rtP->P_794, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T_IWORK_h.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T_IWORK_h.Head] = simTime;
        (*uBuffer)[_rtDW->T_IWORK_h.Head] = _rtB->B_97_338_0;    
  }
             {
    real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK_b.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK_b.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T1_IWORK_f.Head = ((_rtDW->T1_IWORK_f.Head < (_rtDW->T1_IWORK_f.CircularBufSize-1)) ? (_rtDW->T1_IWORK_f.Head+1) : 0);

      if (_rtDW->T1_IWORK_f.Head == _rtDW->T1_IWORK_f.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T1_IWORK_f.CircularBufSize, &_rtDW->T1_IWORK_f.Tail, &_rtDW->T1_IWORK_f.Head, &_rtDW->T1_IWORK_f.Last, simTime - _rtP->P_797, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T1_IWORK_f.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T1_IWORK_f.Head] = simTime;
        (*uBuffer)[_rtDW->T1_IWORK_f.Head] = _rtB->B_97_341_0;    
  }
             {
    real_T **uBuffer = (real_T**)&_rtDW->T_PWORK_b.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T_PWORK_b.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T_IWORK_n.Head = ((_rtDW->T_IWORK_n.Head < (_rtDW->T_IWORK_n.CircularBufSize-1)) ? (_rtDW->T_IWORK_n.Head+1) : 0);

      if (_rtDW->T_IWORK_n.Head == _rtDW->T_IWORK_n.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T_IWORK_n.CircularBufSize, &_rtDW->T_IWORK_n.Tail, &_rtDW->T_IWORK_n.Head, &_rtDW->T_IWORK_n.Last, simTime - _rtP->P_800, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T_IWORK_n.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T_IWORK_n.Head] = simTime;
        (*uBuffer)[_rtDW->T_IWORK_n.Head] = _rtB->B_97_345_0;    
  }
             {
    real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK_j.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK_j.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T1_IWORK_fz.Head = ((_rtDW->T1_IWORK_fz.Head < (_rtDW->T1_IWORK_fz.CircularBufSize-1)) ? (_rtDW->T1_IWORK_fz.Head+1) : 0);

      if (_rtDW->T1_IWORK_fz.Head == _rtDW->T1_IWORK_fz.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T1_IWORK_fz.CircularBufSize, &_rtDW->T1_IWORK_fz.Tail, &_rtDW->T1_IWORK_fz.Head, &_rtDW->T1_IWORK_fz.Last, simTime - _rtP->P_803, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T1_IWORK_fz.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T1_IWORK_fz.Head] = simTime;
        (*uBuffer)[_rtDW->T1_IWORK_fz.Head] = _rtB->B_97_348_0;    
  }
 i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    for (i = 0; i < 6; i++) {
        _rtDW->UnitDelay_DSTATE_ew[i] = _rtB->B_97_370_0[i];
    }
}
i = ssIsSampleHit(S, 3, 0);
if (i != 0) {
    _rtDW->UnitDelay2_DSTATE[0] = _rtB->B_97_508_0[0];
    _rtDW->UnitDelay2_DSTATE[1] = _rtB->B_97_508_0[1];
    _rtDW->UnitDelay2_DSTATE[2] = _rtB->B_97_508_0[2];
}
            {
    real_T **uBuffer = (real_T**)&_rtDW->TransportDelay1_PWORK.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->TransportDelay1_PWORK.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->TransportDelay1_IWORK.Head = ((_rtDW->TransportDelay1_IWORK.Head < (_rtDW->TransportDelay1_IWORK.CircularBufSize-1)) ? (_rtDW->TransportDelay1_IWORK.Head+1) : 0);

      if (_rtDW->TransportDelay1_IWORK.Head == _rtDW->TransportDelay1_IWORK.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->TransportDelay1_IWORK.CircularBufSize, &_rtDW->TransportDelay1_IWORK.Tail, &_rtDW->TransportDelay1_IWORK.Head, &_rtDW->TransportDelay1_IWORK.Last, simTime - _rtP->P_830, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->TransportDelay1_IWORK.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->TransportDelay1_IWORK.Head] = simTime;
        (*uBuffer)[_rtDW->TransportDelay1_IWORK.Head] = _rtB->B_97_405_0;    
  }
 if (_rtDW->GridSupportingasCurrentSourceGridFeeding_MODE) {
    i = ssIsSampleHit(S, 2, 0);
    if (i != 0) {
        _rtDW->DiscreteTimeIntegrator_DSTATE_i += _rtP->P_117 * _rtB->B_27_200_0;
        _rtDW->Integrator_DSTATE_b = _rtP->P_130 * _rtB->B_27_112_0 + _rtB->B_27_113_0;
        _rtDW->Integrator_DSTATE_e = _rtP->P_137 * _rtB->B_27_119_0 + _rtB->B_27_120_0;
        _rtDW->Integrator_DSTATE_b2[0] += _rtP->P_142 * _rtB->B_27_162_0[0];
        _rtDW->Integrator_DSTATE_b2[1] += _rtP->P_142 * _rtB->B_27_162_0[1];
        _rtDW->Integ4_SYSTEM_ENABLE_fj = 0U;
        _rtDW->Integ4_DSTATE_p3 = _rtP->P_164 * _rtB->B_27_175_0 + _rtB->B_27_176_0;
        _rtDW->UnitDelay_DSTATE_oo = _rtB->B_27_206_0;
                    /* Level2 S-Function Block: '<S484>/B_27_17' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 27, 182, SS_CALL_MDL_UPDATE);

         _rtDW->UnitDelay_DSTATE_k2 = _rtB->B_27_175_0;
        if (_rtDW->AutomaticGainControl_MODE_c) {
            i = ssIsSampleHit(S, 2, 0);
            if (i != 0) {
                _rtDW->Integ4_SYSTEM_ENABLE_gq = 0U;
                _rtDW->Integ4_DSTATE_jl = _rtP->P_29 * _rtB->B_23_10_0 + _rtB->B_23_11_0;
                    /* Level2 S-Function Block: '<S479>/B_23_0' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 23, 16, SS_CALL_MDL_UPDATE);

                 _rtDW->UnitDelay_DSTATE_cu = _rtB->B_23_10_0;
                _rtDW->Integ4_SYSTEM_ENABLE_li = 0U;
                _rtDW->Integ4_DSTATE_it = _rtP->P_40 * _rtB->B_23_22_0 + _rtB->B_23_23_0;
                    /* Level2 S-Function Block: '<S481>/B_23_2' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 23, 28, SS_CALL_MDL_UPDATE);

                 _rtDW->UnitDelay_DSTATE_eb = _rtB->B_23_22_0;
            }
        }
        _rtDW->DiscreteTimeIntegrator_DSTATE_oj += _rtP->P_178 * _rtB->B_27_208_0;
        if (_rtDW->DiscreteTimeIntegrator_DSTATE_oj >= _rtP->P_180) {
            _rtDW->DiscreteTimeIntegrator_DSTATE_oj = _rtP->P_180;
        } else {
            if (_rtDW->DiscreteTimeIntegrator_DSTATE_oj <= _rtP->P_181) {
                _rtDW->DiscreteTimeIntegrator_DSTATE_oj = _rtP->P_181;
            }
        }
        _rtDW->UD_DSTATE_i = _rtB->B_27_196_0;
        _rtDW->UnitDelay_DSTATE_b4 = _rtB->B_27_205_0;
              {
          real_T xnew[2];
        xnew[0] = (_rtP->P_191[0])*_rtDW->DiscreteStateSpace_DSTATE_b[0] + (_rtP->P_191[1])*_rtDW->DiscreteStateSpace_DSTATE_b[1];
    xnew[0] += (_rtP->P_192[0])*_rtB->B_27_205_0;

    xnew[1] = (_rtP->P_191[2])*_rtDW->DiscreteStateSpace_DSTATE_b[0] + (_rtP->P_191[3])*_rtDW->DiscreteStateSpace_DSTATE_b[1];
    xnew[1] += (_rtP->P_192[1])*_rtB->B_27_205_0;
    (void) memcpy(&_rtDW->DiscreteStateSpace_DSTATE_b[0], xnew,
sizeof(real_T)*2);
    }

     }
    i = ssIsSampleHit(S, 4, 0);
    if (i != 0) {
        _rtDW->UnitDelay_DSTATE_f3 = _rtB->B_27_234_0;
    }
}
if (_rtDW->Gridforming_MODE) {
    i = ssIsSampleHit(S, 2, 0);
    if (i != 0) {
        _rtDW->DiscreteTimeIntegrator_DSTATE_b += _rtP->P_258 * _rtB->B_36_47_0;
        if (_rtDW->DiscreteTimeIntegrator_DSTATE_b >= _rtP->P_260) {
            _rtDW->DiscreteTimeIntegrator_DSTATE_b = _rtP->P_260;
        } else {
            if (_rtDW->DiscreteTimeIntegrator_DSTATE_b <= _rtP->P_261) {
                _rtDW->DiscreteTimeIntegrator_DSTATE_b = _rtP->P_261;
            }
        }
        _rtDW->DiscreteTimeIntegrator_DSTATE_m += _rtP->P_272 * _rtB->B_36_46_0;
        if (_rtDW->DiscreteTimeIntegrator_DSTATE_m >= _rtP->P_274) {
            _rtDW->DiscreteTimeIntegrator_DSTATE_m = _rtP->P_274;
        } else {
            if (_rtDW->DiscreteTimeIntegrator_DSTATE_m <= _rtP->P_275) {
                _rtDW->DiscreteTimeIntegrator_DSTATE_m = _rtP->P_275;
            }
        }
        _rtDW->DiscreteTimeIntegrator_DSTATE_d += _rtP->P_291 * _rtB->B_36_90_0;
        _rtDW->Integ4_SYSTEM_ENABLE_bq = 0U;
        _rtDW->Integ4_DSTATE_ba = _rtP->P_296 * _rtB->B_36_65_0 + _rtB->B_36_66_0;
        _rtDW->UnitDelay_DSTATE_lw = _rtB->B_36_96_0;
                    /* Level2 S-Function Block: '<S646>/B_36_3' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 36, 72, SS_CALL_MDL_UPDATE);

         _rtDW->UnitDelay_DSTATE_aq = _rtB->B_36_65_0;
        if (_rtDW->AutomaticGainControl_MODE_ab) {
            i = ssIsSampleHit(S, 2, 0);
            if (i != 0) {
                _rtDW->Integ4_SYSTEM_ENABLE_lpf = 0U;
                _rtDW->Integ4_DSTATE_kj = _rtP->P_225 * _rtB->B_30_10_0 + _rtB->B_30_11_0;
                    /* Level2 S-Function Block: '<S641>/B_30_0' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 30, 16, SS_CALL_MDL_UPDATE);

                 _rtDW->UnitDelay_DSTATE_fa = _rtB->B_30_10_0;
                _rtDW->Integ4_SYSTEM_ENABLE_fw = 0U;
                _rtDW->Integ4_DSTATE_d2 = _rtP->P_236 * _rtB->B_30_22_0 + _rtB->B_30_23_0;
                    /* Level2 S-Function Block: '<S643>/B_30_2' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 30, 28, SS_CALL_MDL_UPDATE);

                 _rtDW->UnitDelay_DSTATE_lu = _rtB->B_30_22_0;
            }
        }
        _rtDW->DiscreteTimeIntegrator_DSTATE_ov += _rtP->P_310 * _rtB->B_36_98_0;
        if (_rtDW->DiscreteTimeIntegrator_DSTATE_ov >= _rtP->P_312) {
            _rtDW->DiscreteTimeIntegrator_DSTATE_ov = _rtP->P_312;
        } else {
            if (_rtDW->DiscreteTimeIntegrator_DSTATE_ov <= _rtP->P_313) {
                _rtDW->DiscreteTimeIntegrator_DSTATE_ov = _rtP->P_313;
            }
        }
        _rtDW->UD_DSTATE_e = _rtB->B_36_86_0;
        _rtDW->UnitDelay_DSTATE_n = _rtB->B_36_95_0;
              {
          real_T xnew[2];
        xnew[0] = (_rtP->P_323[0])*_rtDW->DiscreteStateSpace_DSTATE_pm[0] + (_rtP->P_323[1])*_rtDW->DiscreteStateSpace_DSTATE_pm[1];
    xnew[0] += (_rtP->P_324[0])*_rtB->B_36_95_0;

    xnew[1] = (_rtP->P_323[2])*_rtDW->DiscreteStateSpace_DSTATE_pm[0] + (_rtP->P_323[3])*_rtDW->DiscreteStateSpace_DSTATE_pm[1];
    xnew[1] += (_rtP->P_324[1])*_rtB->B_36_95_0;
    (void) memcpy(&_rtDW->DiscreteStateSpace_DSTATE_pm[0], xnew,
sizeof(real_T)*2);
    }

     }
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtDW->DiscreteTimeIntegrator_DSTATE_e += _rtP->P_854 * _rtB->B_97_481_0;
    _rtDW->Integ4_SYSTEM_ENABLE_e = 0U;
    _rtDW->Integ4_DSTATE_p = _rtP->P_859 * _rtB->B_97_456_0 + _rtB->B_97_457_0;
    _rtDW->UnitDelay_DSTATE_i = _rtB->B_97_487_0;
                    /* Level2 S-Function Block: '<S686>/B_82_72' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 463, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_j = _rtB->B_97_456_0;
    if (_rtDW->AutomaticGainControl_MODE_f) {
        i = ssIsSampleHit(S, 2, 0);
        if (i != 0) {
            _rtDW->Integ4_SYSTEM_ENABLE_ht = 0U;
            _rtDW->Integ4_DSTATE_ed = _rtP->P_342 * _rtB->B_44_10_0 + _rtB->B_44_11_0;
                    /* Level2 S-Function Block: '<S681>/B_42_0' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 44, 16, SS_CALL_MDL_UPDATE);

             _rtDW->UnitDelay_DSTATE_g0 = _rtB->B_44_10_0;
            _rtDW->Integ4_SYSTEM_ENABLE_en = 0U;
            _rtDW->Integ4_DSTATE_ok = _rtP->P_353 * _rtB->B_44_22_0 + _rtB->B_44_23_0;
                    /* Level2 S-Function Block: '<S683>/B_42_2' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 44, 28, SS_CALL_MDL_UPDATE);

             _rtDW->UnitDelay_DSTATE_il = _rtB->B_44_22_0;
        }
    }
    _rtDW->DiscreteTimeIntegrator_DSTATE_j += _rtP->P_873 * _rtB->B_97_489_0;
    if (_rtDW->DiscreteTimeIntegrator_DSTATE_j >= _rtP->P_875) {
        _rtDW->DiscreteTimeIntegrator_DSTATE_j = _rtP->P_875;
    } else {
        if (_rtDW->DiscreteTimeIntegrator_DSTATE_j <= _rtP->P_876) {
            _rtDW->DiscreteTimeIntegrator_DSTATE_j = _rtP->P_876;
        }
    }
    _rtDW->UD_DSTATE = _rtB->B_97_477_0;
    _rtDW->UnitDelay_DSTATE_pc = _rtB->B_97_486_0;
              {
          real_T xnew[2];
        xnew[0] = (_rtP->P_886[0])*_rtDW->DiscreteStateSpace_DSTATE[0] + (_rtP->P_886[1])*_rtDW->DiscreteStateSpace_DSTATE[1];
    xnew[0] += (_rtP->P_887[0])*_rtB->B_97_486_0;

    xnew[1] = (_rtP->P_886[2])*_rtDW->DiscreteStateSpace_DSTATE[0] + (_rtP->P_886[3])*_rtDW->DiscreteStateSpace_DSTATE[1];
    xnew[1] += (_rtP->P_887[1])*_rtB->B_97_486_0;
    (void) memcpy(&_rtDW->DiscreteStateSpace_DSTATE[0], xnew,
sizeof(real_T)*2);
    }

     _rtDW->Integ4_SYSTEM_ENABLE_j = 0U;
    _rtDW->Integ4_DSTATE_lr = _rtP->P_910 * _rtB->B_97_586_0 + _rtB->B_97_587_0;
                    /* Level2 S-Function Block: '<S182>/B_82_134' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 589, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_dl = _rtB->B_97_586_0;
    _rtDW->Integ4_SYSTEM_ENABLE_fv = 0U;
    _rtDW->Integ4_DSTATE_g = _rtP->P_920 * _rtB->B_97_596_0 + _rtB->B_97_597_0;
                    /* Level2 S-Function Block: '<S183>/B_82_136' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 599, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_e0 = _rtB->B_97_596_0;
    _rtDW->Integ4_SYSTEM_ENABLE_jr = 0U;
    _rtDW->Integ4_DSTATE_j = _rtP->P_937 * _rtB->B_97_625_0 + _rtB->B_97_626_0;
                    /* Level2 S-Function Block: '<S176>/B_82_138' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 628, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_k = _rtB->B_97_625_0;
    _rtDW->Integ4_SYSTEM_ENABLE_g = 0U;
    _rtDW->Integ4_DSTATE_b = _rtP->P_947 * _rtB->B_97_635_0 + _rtB->B_97_636_0;
                    /* Level2 S-Function Block: '<S177>/B_82_140' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 638, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_g = _rtB->B_97_635_0;
    HoldSine = _rtDW->lastSin_l;
    _rtDW->lastSin_l = _rtDW->lastSin_l * _rtP->P_970 + _rtDW->lastCos_n * _rtP->P_969;
    _rtDW->lastCos_n = _rtDW->lastCos_n * _rtP->P_970 - HoldSine * _rtP->P_969;
    _rtDW->Integ4_SYSTEM_ENABLE_fl = 0U;
    _rtDW->Integ4_DSTATE_a[0] = _rtP->P_973 * _rtB->B_97_665_0[0] + _rtB->B_97_666_0[0];
    _rtDW->Integ4_DSTATE_a[1] = _rtP->P_973 * _rtB->B_97_665_0[1] + _rtB->B_97_666_0[1];
    _rtDW->Integ4_DSTATE_a[2] = _rtP->P_973 * _rtB->B_97_665_0[2] + _rtB->B_97_666_0[2];
                    /* Level2 S-Function Block: '<S778>/B_82_142' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 668, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_g5[0] = _rtB->B_97_665_0[0];
    _rtDW->UnitDelay1_DSTATE_g[0] = _rtB->B_97_675_0[0];
    _rtDW->UnitDelay_DSTATE_g5[1] = _rtB->B_97_665_0[1];
    _rtDW->UnitDelay1_DSTATE_g[1] = _rtB->B_97_675_0[1];
    _rtDW->UnitDelay_DSTATE_g5[2] = _rtB->B_97_665_0[2];
    _rtDW->UnitDelay1_DSTATE_g[2] = _rtB->B_97_675_0[2];
    HoldSine = _rtDW->lastSin_j;
    _rtDW->lastSin_j = _rtDW->lastSin_j * _rtP->P_991 + _rtDW->lastCos_ne * _rtP->P_990;
    _rtDW->lastCos_ne = _rtDW->lastCos_ne * _rtP->P_991 - HoldSine * _rtP->P_990;
    _rtDW->Integ4_SYSTEM_ENABLE_ei = 0U;
    _rtDW->Integ4_DSTATE_e[0] = _rtP->P_994 * _rtB->B_97_677_0[0] + _rtB->B_97_678_0[0];
    _rtDW->Integ4_DSTATE_e[1] = _rtP->P_994 * _rtB->B_97_677_0[1] + _rtB->B_97_678_0[1];
    _rtDW->Integ4_DSTATE_e[2] = _rtP->P_994 * _rtB->B_97_677_0[2] + _rtB->B_97_678_0[2];
                    /* Level2 S-Function Block: '<S776>/B_82_144' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 680, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_b[0] = _rtB->B_97_677_0[0];
    _rtDW->UnitDelay1_DSTATE_j2[0] = _rtB->B_97_687_0[0];
    _rtDW->UnitDelay_DSTATE_b[1] = _rtB->B_97_677_0[1];
    _rtDW->UnitDelay1_DSTATE_j2[1] = _rtB->B_97_687_0[1];
    _rtDW->UnitDelay_DSTATE_b[2] = _rtB->B_97_677_0[2];
    _rtDW->UnitDelay1_DSTATE_j2[2] = _rtB->B_97_687_0[2];
    HoldSine = _rtDW->lastSin_h;
    _rtDW->lastSin_h = _rtDW->lastSin_h * _rtP->P_1016 + _rtDW->lastCos_n3 * _rtP->P_1015;
    _rtDW->lastCos_n3 = _rtDW->lastCos_n3 * _rtP->P_1016 - HoldSine * _rtP->P_1015;
    _rtDW->Integ4_SYSTEM_ENABLE_f0 = 0U;
    _rtDW->Integ4_DSTATE_i[0] = _rtP->P_1019 * _rtB->B_97_695_0[0] + _rtB->B_97_696_0[0];
    _rtDW->Integ4_DSTATE_i[1] = _rtP->P_1019 * _rtB->B_97_695_0[1] + _rtB->B_97_696_0[1];
    _rtDW->Integ4_DSTATE_i[2] = _rtP->P_1019 * _rtB->B_97_695_0[2] + _rtB->B_97_696_0[2];
                    /* Level2 S-Function Block: '<S784>/B_82_146' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 698, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_dlx[0] = _rtB->B_97_695_0[0];
    _rtDW->UnitDelay1_DSTATE_j5[0] = _rtB->B_97_705_0[0];
    _rtDW->UnitDelay_DSTATE_dlx[1] = _rtB->B_97_695_0[1];
    _rtDW->UnitDelay1_DSTATE_j5[1] = _rtB->B_97_705_0[1];
    _rtDW->UnitDelay_DSTATE_dlx[2] = _rtB->B_97_695_0[2];
    _rtDW->UnitDelay1_DSTATE_j5[2] = _rtB->B_97_705_0[2];
    HoldSine = _rtDW->lastSin_n;
    _rtDW->lastSin_n = _rtDW->lastSin_n * _rtP->P_1037 + _rtDW->lastCos_d * _rtP->P_1036;
    _rtDW->lastCos_d = _rtDW->lastCos_d * _rtP->P_1037 - HoldSine * _rtP->P_1036;
    _rtDW->Integ4_SYSTEM_ENABLE_m = 0U;
    _rtDW->Integ4_DSTATE_kd[0] = _rtP->P_1040 * _rtB->B_97_707_0[0] + _rtB->B_97_708_0[0];
    _rtDW->Integ4_DSTATE_kd[1] = _rtP->P_1040 * _rtB->B_97_707_0[1] + _rtB->B_97_708_0[1];
    _rtDW->Integ4_DSTATE_kd[2] = _rtP->P_1040 * _rtB->B_97_707_0[2] + _rtB->B_97_708_0[2];
                    /* Level2 S-Function Block: '<S782>/B_82_148' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 710, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_m[0] = _rtB->B_97_707_0[0];
    _rtDW->UnitDelay1_DSTATE_a[0] = _rtB->B_97_717_0[0];
    _rtDW->UnitDelay_DSTATE_m[1] = _rtB->B_97_707_0[1];
    _rtDW->UnitDelay1_DSTATE_a[1] = _rtB->B_97_717_0[1];
    _rtDW->UnitDelay_DSTATE_m[2] = _rtB->B_97_707_0[2];
    _rtDW->UnitDelay1_DSTATE_a[2] = _rtB->B_97_717_0[2];
}
            {
    real_T **uBuffer = (real_T**)&_rtDW->T_PWORK_bx.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T_PWORK_bx.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T_IWORK_a.Head = ((_rtDW->T_IWORK_a.Head < (_rtDW->T_IWORK_a.CircularBufSize-1)) ? (_rtDW->T_IWORK_a.Head+1) : 0);

      if (_rtDW->T_IWORK_a.Head == _rtDW->T_IWORK_a.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T_IWORK_a.CircularBufSize, &_rtDW->T_IWORK_a.Tail, &_rtDW->T_IWORK_a.Head, &_rtDW->T_IWORK_a.Last, simTime - _rtP->P_1059, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T_IWORK_a.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T_IWORK_a.Head] = simTime;
        (*uBuffer)[_rtDW->T_IWORK_a.Head] = _rtB->B_97_729_0;    
  }
             {
    real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK_o.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK_o.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T1_IWORK_n.Head = ((_rtDW->T1_IWORK_n.Head < (_rtDW->T1_IWORK_n.CircularBufSize-1)) ? (_rtDW->T1_IWORK_n.Head+1) : 0);

      if (_rtDW->T1_IWORK_n.Head == _rtDW->T1_IWORK_n.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T1_IWORK_n.CircularBufSize, &_rtDW->T1_IWORK_n.Tail, &_rtDW->T1_IWORK_n.Head, &_rtDW->T1_IWORK_n.Last, simTime - _rtP->P_1062, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T1_IWORK_n.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T1_IWORK_n.Head] = simTime;
        (*uBuffer)[_rtDW->T1_IWORK_n.Head] = _rtB->B_97_732_0;    
  }
             {
    real_T **uBuffer = (real_T**)&_rtDW->T_PWORK_le.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T_PWORK_le.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T_IWORK_c.Head = ((_rtDW->T_IWORK_c.Head < (_rtDW->T_IWORK_c.CircularBufSize-1)) ? (_rtDW->T_IWORK_c.Head+1) : 0);

      if (_rtDW->T_IWORK_c.Head == _rtDW->T_IWORK_c.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T_IWORK_c.CircularBufSize, &_rtDW->T_IWORK_c.Tail, &_rtDW->T_IWORK_c.Head, &_rtDW->T_IWORK_c.Last, simTime - _rtP->P_1065, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T_IWORK_c.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T_IWORK_c.Head] = simTime;
        (*uBuffer)[_rtDW->T_IWORK_c.Head] = _rtB->B_97_736_0;    
  }
             {
    real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK_c.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK_c.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T1_IWORK_c.Head = ((_rtDW->T1_IWORK_c.Head < (_rtDW->T1_IWORK_c.CircularBufSize-1)) ? (_rtDW->T1_IWORK_c.Head+1) : 0);

      if (_rtDW->T1_IWORK_c.Head == _rtDW->T1_IWORK_c.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T1_IWORK_c.CircularBufSize, &_rtDW->T1_IWORK_c.Tail, &_rtDW->T1_IWORK_c.Head, &_rtDW->T1_IWORK_c.Last, simTime - _rtP->P_1068, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T1_IWORK_c.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T1_IWORK_c.Head] = simTime;
        (*uBuffer)[_rtDW->T1_IWORK_c.Head] = _rtB->B_97_739_0;    
  }
             {
    real_T **uBuffer = (real_T**)&_rtDW->T_PWORK_i.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T_PWORK_i.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T_IWORK_o.Head = ((_rtDW->T_IWORK_o.Head < (_rtDW->T_IWORK_o.CircularBufSize-1)) ? (_rtDW->T_IWORK_o.Head+1) : 0);

      if (_rtDW->T_IWORK_o.Head == _rtDW->T_IWORK_o.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T_IWORK_o.CircularBufSize, &_rtDW->T_IWORK_o.Tail, &_rtDW->T_IWORK_o.Head, &_rtDW->T_IWORK_o.Last, simTime - _rtP->P_1071, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T_IWORK_o.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T_IWORK_o.Head] = simTime;
        (*uBuffer)[_rtDW->T_IWORK_o.Head] = _rtB->B_97_743_0;    
  }
             {
    real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK_h.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK_h.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T1_IWORK_j.Head = ((_rtDW->T1_IWORK_j.Head < (_rtDW->T1_IWORK_j.CircularBufSize-1)) ? (_rtDW->T1_IWORK_j.Head+1) : 0);

      if (_rtDW->T1_IWORK_j.Head == _rtDW->T1_IWORK_j.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T1_IWORK_j.CircularBufSize, &_rtDW->T1_IWORK_j.Tail, &_rtDW->T1_IWORK_j.Head, &_rtDW->T1_IWORK_j.Last, simTime - _rtP->P_1074, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T1_IWORK_j.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T1_IWORK_j.Head] = simTime;
        (*uBuffer)[_rtDW->T1_IWORK_j.Head] = _rtB->B_97_746_0;    
  }
 i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    for (i = 0; i < 6; i++) {
        _rtDW->UnitDelay_DSTATE_c[i] = _rtB->B_97_768_0[i];
    }
}
i = ssIsSampleHit(S, 3, 0);
if (i != 0) {
    _rtDW->UnitDelay4_DSTATE[0] = _rtB->B_97_1037_0[0];
    _rtDW->UnitDelay4_DSTATE[1] = _rtB->B_97_1037_0[1];
    _rtDW->UnitDelay4_DSTATE[2] = _rtB->B_97_1037_0[2];
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtDW->DiscreteTimeIntegrator_DSTATE_p += _rtP->P_1119 * _rtB->B_97_867_0;
    _rtDW->Integ4_SYSTEM_ENABLE_i = 0U;
    _rtDW->Integ4_DSTATE_ay = _rtP->P_1124 * _rtB->B_97_842_0 + _rtB->B_97_843_0;
    _rtDW->UnitDelay_DSTATE_l = _rtB->B_97_873_0;
                    /* Level2 S-Function Block: '<S770>/B_82_169' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 849, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_dw = _rtB->B_97_842_0;
    if (_rtDW->AutomaticGainControl_MODE_p) {
        i = ssIsSampleHit(S, 2, 0);
        if (i != 0) {
            _rtDW->Integ4_SYSTEM_ENABLE_k5 = 0U;
            _rtDW->Integ4_DSTATE_o = _rtP->P_422 * _rtB->B_66_10_0 + _rtB->B_66_11_0;
                    /* Level2 S-Function Block: '<S765>/B_58_0' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 66, 16, SS_CALL_MDL_UPDATE);

             _rtDW->UnitDelay_DSTATE_a = _rtB->B_66_10_0;
            _rtDW->Integ4_SYSTEM_ENABLE_ks = 0U;
            _rtDW->Integ4_DSTATE_g1 = _rtP->P_433 * _rtB->B_66_22_0 + _rtB->B_66_23_0;
                    /* Level2 S-Function Block: '<S767>/B_58_2' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 66, 28, SS_CALL_MDL_UPDATE);

             _rtDW->UnitDelay_DSTATE_b0 = _rtB->B_66_22_0;
        }
    }
    _rtDW->DiscreteTimeIntegrator_DSTATE_pw += _rtP->P_1138 * _rtB->B_97_875_0;
    if (_rtDW->DiscreteTimeIntegrator_DSTATE_pw >= _rtP->P_1140) {
        _rtDW->DiscreteTimeIntegrator_DSTATE_pw = _rtP->P_1140;
    } else {
        if (_rtDW->DiscreteTimeIntegrator_DSTATE_pw <= _rtP->P_1141) {
            _rtDW->DiscreteTimeIntegrator_DSTATE_pw = _rtP->P_1141;
        }
    }
    _rtDW->UD_DSTATE_b = _rtB->B_97_863_0;
    _rtDW->UnitDelay_DSTATE_kd = _rtB->B_97_872_0;
              {
          real_T xnew[2];
        xnew[0] = (_rtP->P_1151[0])*_rtDW->DiscreteStateSpace_DSTATE_i[0] + (_rtP->P_1151[1])*_rtDW->DiscreteStateSpace_DSTATE_i[1];
    xnew[0] += (_rtP->P_1152[0])*_rtB->B_97_872_0;

    xnew[1] = (_rtP->P_1151[2])*_rtDW->DiscreteStateSpace_DSTATE_i[0] + (_rtP->P_1151[3])*_rtDW->DiscreteStateSpace_DSTATE_i[1];
    xnew[1] += (_rtP->P_1152[1])*_rtB->B_97_872_0;
    (void) memcpy(&_rtDW->DiscreteStateSpace_DSTATE_i[0], xnew,
sizeof(real_T)*2);
    }

     _rtDW->DiscreteTimeIntegrator_DSTATE_o += _rtP->P_1220 * _rtB->B_97_1085_0;
    _rtDW->Integrator_DSTATE = _rtP->P_1233 * _rtB->B_97_998_0 + _rtB->B_97_999_0;
    _rtDW->Integrator_DSTATE_l = _rtP->P_1240 * _rtB->B_97_1005_0 + _rtB->B_97_1006_0;
    _rtDW->Integrator_DSTATE_n[0] += _rtP->P_1245 * _rtB->B_97_1048_0[0];
    _rtDW->Integrator_DSTATE_n[1] += _rtP->P_1245 * _rtB->B_97_1048_0[1];
    _rtDW->Integ4_SYSTEM_ENABLE_hq = 0U;
    _rtDW->Integ4_DSTATE_bi = _rtP->P_1267 * _rtB->B_97_1060_0 + _rtB->B_97_1061_0;
    _rtDW->UnitDelay_DSTATE_dd = _rtB->B_97_1091_0;
                    /* Level2 S-Function Block: '<S1155>/B_82_190' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 1067, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_mf = _rtB->B_97_1060_0;
    if (_rtDW->AutomaticGainControl_MODE) {
        i = ssIsSampleHit(S, 2, 0);
        if (i != 0) {
            _rtDW->Integ4_SYSTEM_ENABLE_f0v = 0U;
            _rtDW->Integ4_DSTATE_pa = _rtP->P_508 * _rtB->B_93_10_0 + _rtB->B_93_11_0;
                    /* Level2 S-Function Block: '<S1150>/B_78_0' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 93, 16, SS_CALL_MDL_UPDATE);

             _rtDW->UnitDelay_DSTATE_eer = _rtB->B_93_10_0;
            _rtDW->Integ4_SYSTEM_ENABLE_ki = 0U;
            _rtDW->Integ4_DSTATE_jc = _rtP->P_519 * _rtB->B_93_22_0 + _rtB->B_93_23_0;
                    /* Level2 S-Function Block: '<S1152>/B_78_2' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 93, 28, SS_CALL_MDL_UPDATE);

             _rtDW->UnitDelay_DSTATE_fv = _rtB->B_93_22_0;
        }
    }
    _rtDW->DiscreteTimeIntegrator_DSTATE_jp += _rtP->P_1281 * _rtB->B_97_1093_0;
    if (_rtDW->DiscreteTimeIntegrator_DSTATE_jp >= _rtP->P_1283) {
        _rtDW->DiscreteTimeIntegrator_DSTATE_jp = _rtP->P_1283;
    } else {
        if (_rtDW->DiscreteTimeIntegrator_DSTATE_jp <= _rtP->P_1284) {
            _rtDW->DiscreteTimeIntegrator_DSTATE_jp = _rtP->P_1284;
        }
    }
    _rtDW->UD_DSTATE_a = _rtB->B_97_1081_0;
    _rtDW->UnitDelay_DSTATE_kv = _rtB->B_97_1090_0;
              {
          real_T xnew[2];
        xnew[0] = (_rtP->P_1294[0])*_rtDW->DiscreteStateSpace_DSTATE_p[0] + (_rtP->P_1294[1])*_rtDW->DiscreteStateSpace_DSTATE_p[1];
    xnew[0] += (_rtP->P_1295[0])*_rtB->B_97_1090_0;

    xnew[1] = (_rtP->P_1294[2])*_rtDW->DiscreteStateSpace_DSTATE_p[0] + (_rtP->P_1294[3])*_rtDW->DiscreteStateSpace_DSTATE_p[1];
    xnew[1] += (_rtP->P_1295[1])*_rtB->B_97_1090_0;
    (void) memcpy(&_rtDW->DiscreteStateSpace_DSTATE_p[0], xnew,
sizeof(real_T)*2);
    }

 }
i = ssIsSampleHit(S, 4, 0);
if (i != 0) {
    _rtDW->UnitDelay_DSTATE_ef = _rtB->B_97_1119_0;
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
              {
          real_T xnew[1];
        xnew[0] = _rtP->P_1339*_rtDW->DiscreteStateSpace_DSTATE_n;
    xnew[0] += _rtP->P_1340*_rtB->B_97_1196_0;
    (void) memcpy(&_rtDW->DiscreteStateSpace_DSTATE_n, xnew,
sizeof(real_T)*1);
    }

     _rtDW->UnitDelay1_DSTATE_j5g = _rtB->B_97_1226_0;
              {
          real_T xnew[1];
        xnew[0] = _rtP->P_1351*_rtDW->DiscreteStateSpace_DSTATE_m;
    xnew[0] += _rtP->P_1352*_rtB->B_97_1201_0;
    (void) memcpy(&_rtDW->DiscreteStateSpace_DSTATE_m, xnew,
sizeof(real_T)*1);
    }

               {
          real_T xnew[1];
        xnew[0] = (_rtP->P_1362)*_rtDW->DiscreteStateSpace_DSTATE_k;
    xnew[0] += _rtP->P_1363*_rtB->B_97_1223_0;
    (void) memcpy(&_rtDW->DiscreteStateSpace_DSTATE_k, xnew,
sizeof(real_T)*1);
    }

               {
          real_T xnew[1];
        xnew[0] = _rtP->P_1368*_rtDW->DiscreteStateSpace_DSTATE_h;
    xnew[0] += _rtP->P_1369*_rtB->B_97_1225_0;
    (void) memcpy(&_rtDW->DiscreteStateSpace_DSTATE_h, xnew,
sizeof(real_T)*1);
    }

     _rtDW->UnitDelay2_DSTATE_n = _rtB->B_97_1224_0;
}
            {
    real_T **uBuffer = (real_T**)&_rtDW->ENGINETd_PWORK.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->ENGINETd_PWORK.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->ENGINETd_IWORK.Head = ((_rtDW->ENGINETd_IWORK.Head < (_rtDW->ENGINETd_IWORK.CircularBufSize-1)) ? (_rtDW->ENGINETd_IWORK.Head+1) : 0);

      if (_rtDW->ENGINETd_IWORK.Head == _rtDW->ENGINETd_IWORK.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->ENGINETd_IWORK.CircularBufSize, &_rtDW->ENGINETd_IWORK.Tail, &_rtDW->ENGINETd_IWORK.Head, &_rtDW->ENGINETd_IWORK.Last, simTime - _rtP->P_1376, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->ENGINETd_IWORK.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->ENGINETd_IWORK.Head] = simTime;
        (*uBuffer)[_rtDW->ENGINETd_IWORK.Head] = _rtB->B_97_1232_0;    
  }
 /* 0: INTG_NORMAL     1: INTG_LEAVING_UPPER_SAT  2: INTG_LEAVING_LOWER_SAT */
/* 3: INTG_UPPER_SAT  4: INTG_LOWER_SAT */
if (_rtX->Integrator_CSTATE_g1 == _rtP->P_1380) {
    switch (_rtDW->Integrator_MODE) {
      case 3:
        if (_rtB->B_97_1248_0 < 0.0) {
            ssSetBlockStateForSolverChangedAtMajorStep(S);
            _rtDW->Integrator_MODE = 1;
        }
        break;
      case 1:
        if (_rtB->B_97_1248_0 >= 0.0) {
            _rtDW->Integrator_MODE = 3;
            ssSetBlockStateForSolverChangedAtMajorStep(S);
        }
        break;
      default:
        ssSetBlockStateForSolverChangedAtMajorStep(S);
        if (_rtB->B_97_1248_0 < 0.0) {
            _rtDW->Integrator_MODE = 1;
        } else {
            _rtDW->Integrator_MODE = 3;
        }
        break;
    }
} else if (_rtX->Integrator_CSTATE_g1 == _rtP->P_1381) {
    switch (_rtDW->Integrator_MODE) {
      case 4:
        if (_rtB->B_97_1248_0 > 0.0) {
            ssSetBlockStateForSolverChangedAtMajorStep(S);
            _rtDW->Integrator_MODE = 2;
        }
        break;
      case 2:
        if (_rtB->B_97_1248_0 <= 0.0) {
            _rtDW->Integrator_MODE = 4;
            ssSetBlockStateForSolverChangedAtMajorStep(S);
        }
        break;
      default:
        ssSetBlockStateForSolverChangedAtMajorStep(S);
        if (_rtB->B_97_1248_0 > 0.0) {
            _rtDW->Integrator_MODE = 2;
        } else {
            _rtDW->Integrator_MODE = 4;
        }
        break;
    }
} else {
    _rtDW->Integrator_MODE = 0;
}
i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtDW->DiscreteTimeIntegrator2_DSTATE += _rtP->P_1385 * _rtB->B_97_1259_0;
    HoldSine = _rtDW->lastSin_b;
    _rtDW->lastSin_b = _rtDW->lastSin_b * _rtP->P_1410 + _rtDW->lastCos_p * _rtP->P_1409;
    _rtDW->lastCos_p = _rtDW->lastCos_p * _rtP->P_1410 - HoldSine * _rtP->P_1409;
    _rtDW->Integ4_SYSTEM_ENABLE_c2 = 0U;
    _rtDW->Integ4_DSTATE_dr[0] = _rtP->P_1413 * _rtB->B_97_1270_0[0] + _rtB->B_97_1271_0[0];
    _rtDW->Integ4_DSTATE_dr[1] = _rtP->P_1413 * _rtB->B_97_1270_0[1] + _rtB->B_97_1271_0[1];
    _rtDW->Integ4_DSTATE_dr[2] = _rtP->P_1413 * _rtB->B_97_1270_0[2] + _rtB->B_97_1271_0[2];
                    /* Level2 S-Function Block: '<S736>/B_82_239' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 1273, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_is[0] = _rtB->B_97_1270_0[0];
    _rtDW->UnitDelay1_DSTATE_cp[0] = _rtB->B_97_1280_0[0];
    _rtDW->UnitDelay_DSTATE_is[1] = _rtB->B_97_1270_0[1];
    _rtDW->UnitDelay1_DSTATE_cp[1] = _rtB->B_97_1280_0[1];
    _rtDW->UnitDelay_DSTATE_is[2] = _rtB->B_97_1270_0[2];
    _rtDW->UnitDelay1_DSTATE_cp[2] = _rtB->B_97_1280_0[2];
    HoldSine = _rtDW->lastSin_jg;
    _rtDW->lastSin_jg = _rtDW->lastSin_jg * _rtP->P_1431 + _rtDW->lastCos_fb * _rtP->P_1430;
    _rtDW->lastCos_fb = _rtDW->lastCos_fb * _rtP->P_1431 - HoldSine * _rtP->P_1430;
    _rtDW->Integ4_SYSTEM_ENABLE_jh = 0U;
    _rtDW->Integ4_DSTATE_g5[0] = _rtP->P_1434 * _rtB->B_97_1282_0[0] + _rtB->B_97_1283_0[0];
    _rtDW->Integ4_DSTATE_g5[1] = _rtP->P_1434 * _rtB->B_97_1282_0[1] + _rtB->B_97_1283_0[1];
    _rtDW->Integ4_DSTATE_g5[2] = _rtP->P_1434 * _rtB->B_97_1282_0[2] + _rtB->B_97_1283_0[2];
                    /* Level2 S-Function Block: '<S734>/B_82_241' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 1285, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_ox[0] = _rtB->B_97_1282_0[0];
    _rtDW->UnitDelay1_DSTATE_f[0] = _rtB->B_97_1292_0[0];
    _rtDW->UnitDelay_DSTATE_ox[1] = _rtB->B_97_1282_0[1];
    _rtDW->UnitDelay1_DSTATE_f[1] = _rtB->B_97_1292_0[1];
    _rtDW->UnitDelay_DSTATE_ox[2] = _rtB->B_97_1282_0[2];
    _rtDW->UnitDelay1_DSTATE_f[2] = _rtB->B_97_1292_0[2];
    HoldSine = _rtDW->lastSin_g;
    _rtDW->lastSin_g = _rtDW->lastSin_g * _rtP->P_1456 + _rtDW->lastCos_b * _rtP->P_1455;
    _rtDW->lastCos_b = _rtDW->lastCos_b * _rtP->P_1456 - HoldSine * _rtP->P_1455;
    _rtDW->Integ4_SYSTEM_ENABLE_n = 0U;
    _rtDW->Integ4_DSTATE_h3[0] = _rtP->P_1459 * _rtB->B_97_1300_0[0] + _rtB->B_97_1301_0[0];
    _rtDW->Integ4_DSTATE_h3[1] = _rtP->P_1459 * _rtB->B_97_1300_0[1] + _rtB->B_97_1301_0[1];
    _rtDW->Integ4_DSTATE_h3[2] = _rtP->P_1459 * _rtB->B_97_1300_0[2] + _rtB->B_97_1301_0[2];
                    /* Level2 S-Function Block: '<S742>/B_82_243' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 1303, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_bp[0] = _rtB->B_97_1300_0[0];
    _rtDW->UnitDelay1_DSTATE_gq[0] = _rtB->B_97_1310_0[0];
    _rtDW->UnitDelay_DSTATE_bp[1] = _rtB->B_97_1300_0[1];
    _rtDW->UnitDelay1_DSTATE_gq[1] = _rtB->B_97_1310_0[1];
    _rtDW->UnitDelay_DSTATE_bp[2] = _rtB->B_97_1300_0[2];
    _rtDW->UnitDelay1_DSTATE_gq[2] = _rtB->B_97_1310_0[2];
    HoldSine = _rtDW->lastSin_i;
    _rtDW->lastSin_i = _rtDW->lastSin_i * _rtP->P_1477 + _rtDW->lastCos_e * _rtP->P_1476;
    _rtDW->lastCos_e = _rtDW->lastCos_e * _rtP->P_1477 - HoldSine * _rtP->P_1476;
    _rtDW->Integ4_SYSTEM_ENABLE_cn = 0U;
    _rtDW->Integ4_DSTATE_bj[0] = _rtP->P_1480 * _rtB->B_97_1312_0[0] + _rtB->B_97_1313_0[0];
    _rtDW->Integ4_DSTATE_bj[1] = _rtP->P_1480 * _rtB->B_97_1312_0[1] + _rtB->B_97_1313_0[1];
    _rtDW->Integ4_DSTATE_bj[2] = _rtP->P_1480 * _rtB->B_97_1312_0[2] + _rtB->B_97_1313_0[2];
                    /* Level2 S-Function Block: '<S740>/B_82_245' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 1315, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_po[0] = _rtB->B_97_1312_0[0];
    _rtDW->UnitDelay1_DSTATE_gu[0] = _rtB->B_97_1322_0[0];
    _rtDW->UnitDelay_DSTATE_po[1] = _rtB->B_97_1312_0[1];
    _rtDW->UnitDelay1_DSTATE_gu[1] = _rtB->B_97_1322_0[1];
    _rtDW->UnitDelay_DSTATE_po[2] = _rtB->B_97_1312_0[2];
    _rtDW->UnitDelay1_DSTATE_gu[2] = _rtB->B_97_1322_0[2];
}
            {
    real_T **uBuffer = (real_T**)&_rtDW->T_PWORK_g.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T_PWORK_g.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T_IWORK_b.Head = ((_rtDW->T_IWORK_b.Head < (_rtDW->T_IWORK_b.CircularBufSize-1)) ? (_rtDW->T_IWORK_b.Head+1) : 0);

      if (_rtDW->T_IWORK_b.Head == _rtDW->T_IWORK_b.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T_IWORK_b.CircularBufSize, &_rtDW->T_IWORK_b.Tail, &_rtDW->T_IWORK_b.Head, &_rtDW->T_IWORK_b.Last, simTime - _rtP->P_1499, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T_IWORK_b.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T_IWORK_b.Head] = simTime;
        (*uBuffer)[_rtDW->T_IWORK_b.Head] = _rtB->B_97_1334_0;    
  }
             {
    real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK_bd.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK_bd.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T1_IWORK_fg.Head = ((_rtDW->T1_IWORK_fg.Head < (_rtDW->T1_IWORK_fg.CircularBufSize-1)) ? (_rtDW->T1_IWORK_fg.Head+1) : 0);

      if (_rtDW->T1_IWORK_fg.Head == _rtDW->T1_IWORK_fg.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T1_IWORK_fg.CircularBufSize, &_rtDW->T1_IWORK_fg.Tail, &_rtDW->T1_IWORK_fg.Head, &_rtDW->T1_IWORK_fg.Last, simTime - _rtP->P_1502, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T1_IWORK_fg.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T1_IWORK_fg.Head] = simTime;
        (*uBuffer)[_rtDW->T1_IWORK_fg.Head] = _rtB->B_97_1337_0;    
  }
             {
    real_T **uBuffer = (real_T**)&_rtDW->T_PWORK_d.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T_PWORK_d.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T_IWORK_m.Head = ((_rtDW->T_IWORK_m.Head < (_rtDW->T_IWORK_m.CircularBufSize-1)) ? (_rtDW->T_IWORK_m.Head+1) : 0);

      if (_rtDW->T_IWORK_m.Head == _rtDW->T_IWORK_m.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T_IWORK_m.CircularBufSize, &_rtDW->T_IWORK_m.Tail, &_rtDW->T_IWORK_m.Head, &_rtDW->T_IWORK_m.Last, simTime - _rtP->P_1505, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T_IWORK_m.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T_IWORK_m.Head] = simTime;
        (*uBuffer)[_rtDW->T_IWORK_m.Head] = _rtB->B_97_1341_0;    
  }
             {
    real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK_bk.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK_bk.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T1_IWORK_m.Head = ((_rtDW->T1_IWORK_m.Head < (_rtDW->T1_IWORK_m.CircularBufSize-1)) ? (_rtDW->T1_IWORK_m.Head+1) : 0);

      if (_rtDW->T1_IWORK_m.Head == _rtDW->T1_IWORK_m.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T1_IWORK_m.CircularBufSize, &_rtDW->T1_IWORK_m.Tail, &_rtDW->T1_IWORK_m.Head, &_rtDW->T1_IWORK_m.Last, simTime - _rtP->P_1508, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T1_IWORK_m.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T1_IWORK_m.Head] = simTime;
        (*uBuffer)[_rtDW->T1_IWORK_m.Head] = _rtB->B_97_1344_0;    
  }
             {
    real_T **uBuffer = (real_T**)&_rtDW->T_PWORK_l4.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T_PWORK_l4.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T_IWORK_c1.Head = ((_rtDW->T_IWORK_c1.Head < (_rtDW->T_IWORK_c1.CircularBufSize-1)) ? (_rtDW->T_IWORK_c1.Head+1) : 0);

      if (_rtDW->T_IWORK_c1.Head == _rtDW->T_IWORK_c1.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T_IWORK_c1.CircularBufSize, &_rtDW->T_IWORK_c1.Tail, &_rtDW->T_IWORK_c1.Head, &_rtDW->T_IWORK_c1.Last, simTime - _rtP->P_1511, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T_IWORK_c1.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T_IWORK_c1.Head] = simTime;
        (*uBuffer)[_rtDW->T_IWORK_c1.Head] = _rtB->B_97_1348_0;    
  }
             {
    real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK_p.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK_p.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T1_IWORK_c4.Head = ((_rtDW->T1_IWORK_c4.Head < (_rtDW->T1_IWORK_c4.CircularBufSize-1)) ? (_rtDW->T1_IWORK_c4.Head+1) : 0);

      if (_rtDW->T1_IWORK_c4.Head == _rtDW->T1_IWORK_c4.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T1_IWORK_c4.CircularBufSize, &_rtDW->T1_IWORK_c4.Tail, &_rtDW->T1_IWORK_c4.Head, &_rtDW->T1_IWORK_c4.Last, simTime - _rtP->P_1514, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T1_IWORK_c4.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T1_IWORK_c4.Head] = simTime;
        (*uBuffer)[_rtDW->T1_IWORK_c4.Head] = _rtB->B_97_1351_0;    
  }
 i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtDW->DiscreteTimeIntegrator_DSTATE_fi += _rtP->P_1538 * _rtB->B_97_1430_0;
    _rtDW->Integ4_SYSTEM_ENABLE_kv = 0U;
    _rtDW->Integ4_DSTATE_f2 = _rtP->P_1543 * _rtB->B_97_1405_0 + _rtB->B_97_1406_0;
    _rtDW->UnitDelay_DSTATE_l3 = _rtB->B_97_1436_0;
                    /* Level2 S-Function Block: '<S728>/B_82_260' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 1412, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_gl = _rtB->B_97_1405_0;
    if (_rtDW->AutomaticGainControl_MODE_a) {
        i = ssIsSampleHit(S, 2, 0);
        if (i != 0) {
            _rtDW->Integ4_SYSTEM_ENABLE_fs = 0U;
            _rtDW->Integ4_DSTATE_i0 = _rtP->P_382 * _rtB->B_55_10_0 + _rtB->B_55_11_0;
                    /* Level2 S-Function Block: '<S723>/B_50_0' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 55, 16, SS_CALL_MDL_UPDATE);

             _rtDW->UnitDelay_DSTATE_pv = _rtB->B_55_10_0;
            _rtDW->Integ4_SYSTEM_ENABLE_mo = 0U;
            _rtDW->Integ4_DSTATE_c = _rtP->P_393 * _rtB->B_55_22_0 + _rtB->B_55_23_0;
                    /* Level2 S-Function Block: '<S725>/B_50_2' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 55, 28, SS_CALL_MDL_UPDATE);

             _rtDW->UnitDelay_DSTATE_as = _rtB->B_55_22_0;
        }
    }
    _rtDW->DiscreteTimeIntegrator_DSTATE_a += _rtP->P_1557 * _rtB->B_97_1438_0;
    if (_rtDW->DiscreteTimeIntegrator_DSTATE_a >= _rtP->P_1559) {
        _rtDW->DiscreteTimeIntegrator_DSTATE_a = _rtP->P_1559;
    } else {
        if (_rtDW->DiscreteTimeIntegrator_DSTATE_a <= _rtP->P_1560) {
            _rtDW->DiscreteTimeIntegrator_DSTATE_a = _rtP->P_1560;
        }
    }
    _rtDW->UD_DSTATE_aa = _rtB->B_97_1426_0;
    _rtDW->UnitDelay_DSTATE_oi = _rtB->B_97_1435_0;
              {
          real_T xnew[2];
        xnew[0] = (_rtP->P_1570[0])*_rtDW->DiscreteStateSpace_DSTATE_a[0] + (_rtP->P_1570[1])*_rtDW->DiscreteStateSpace_DSTATE_a[1];
    xnew[0] += (_rtP->P_1571[0])*_rtB->B_97_1435_0;

    xnew[1] = (_rtP->P_1570[2])*_rtDW->DiscreteStateSpace_DSTATE_a[0] + (_rtP->P_1570[3])*_rtDW->DiscreteStateSpace_DSTATE_a[1];
    xnew[1] += (_rtP->P_1571[1])*_rtB->B_97_1435_0;
    (void) memcpy(&_rtDW->DiscreteStateSpace_DSTATE_a[0], xnew,
sizeof(real_T)*2);
    }

     HoldSine = _rtDW->lastSin_im;
    _rtDW->lastSin_im = _rtDW->lastSin_im * _rtP->P_1585 + _rtDW->lastCos_k * _rtP->P_1584;
    _rtDW->lastCos_k = _rtDW->lastCos_k * _rtP->P_1585 - HoldSine * _rtP->P_1584;
    _rtDW->Integ4_SYSTEM_ENABLE_pe = 0U;
    _rtDW->Integ4_DSTATE_gh[0] = _rtP->P_1588 * _rtB->B_97_1481_0[0] + _rtB->B_97_1482_0[0];
    _rtDW->Integ4_DSTATE_gh[1] = _rtP->P_1588 * _rtB->B_97_1481_0[1] + _rtB->B_97_1482_0[1];
    _rtDW->Integ4_DSTATE_gh[2] = _rtP->P_1588 * _rtB->B_97_1481_0[2] + _rtB->B_97_1482_0[2];
                    /* Level2 S-Function Block: '<S820>/B_82_298' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 1484, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_bw[0] = _rtB->B_97_1481_0[0];
    _rtDW->UnitDelay1_DSTATE_hp[0] = _rtB->B_97_1491_0[0];
    _rtDW->UnitDelay_DSTATE_bw[1] = _rtB->B_97_1481_0[1];
    _rtDW->UnitDelay1_DSTATE_hp[1] = _rtB->B_97_1491_0[1];
    _rtDW->UnitDelay_DSTATE_bw[2] = _rtB->B_97_1481_0[2];
    _rtDW->UnitDelay1_DSTATE_hp[2] = _rtB->B_97_1491_0[2];
    HoldSine = _rtDW->lastSin_bt;
    _rtDW->lastSin_bt = _rtDW->lastSin_bt * _rtP->P_1606 + _rtDW->lastCos_m * _rtP->P_1605;
    _rtDW->lastCos_m = _rtDW->lastCos_m * _rtP->P_1606 - HoldSine * _rtP->P_1605;
    _rtDW->Integ4_SYSTEM_ENABLE_ix = 0U;
    _rtDW->Integ4_DSTATE_n4[0] = _rtP->P_1609 * _rtB->B_97_1493_0[0] + _rtB->B_97_1494_0[0];
    _rtDW->Integ4_DSTATE_n4[1] = _rtP->P_1609 * _rtB->B_97_1493_0[1] + _rtB->B_97_1494_0[1];
    _rtDW->Integ4_DSTATE_n4[2] = _rtP->P_1609 * _rtB->B_97_1493_0[2] + _rtB->B_97_1494_0[2];
                    /* Level2 S-Function Block: '<S818>/B_82_300' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 1496, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_kh[0] = _rtB->B_97_1493_0[0];
    _rtDW->UnitDelay1_DSTATE_fp[0] = _rtB->B_97_1503_0[0];
    _rtDW->UnitDelay_DSTATE_kh[1] = _rtB->B_97_1493_0[1];
    _rtDW->UnitDelay1_DSTATE_fp[1] = _rtB->B_97_1503_0[1];
    _rtDW->UnitDelay_DSTATE_kh[2] = _rtB->B_97_1493_0[2];
    _rtDW->UnitDelay1_DSTATE_fp[2] = _rtB->B_97_1503_0[2];
    HoldSine = _rtDW->lastSin_p;
    _rtDW->lastSin_p = _rtDW->lastSin_p * _rtP->P_1631 + _rtDW->lastCos_fk * _rtP->P_1630;
    _rtDW->lastCos_fk = _rtDW->lastCos_fk * _rtP->P_1631 - HoldSine * _rtP->P_1630;
    _rtDW->Integ4_SYSTEM_ENABLE_lx = 0U;
    _rtDW->Integ4_DSTATE_pq[0] = _rtP->P_1634 * _rtB->B_97_1511_0[0] + _rtB->B_97_1512_0[0];
    _rtDW->Integ4_DSTATE_pq[1] = _rtP->P_1634 * _rtB->B_97_1511_0[1] + _rtB->B_97_1512_0[1];
    _rtDW->Integ4_DSTATE_pq[2] = _rtP->P_1634 * _rtB->B_97_1511_0[2] + _rtB->B_97_1512_0[2];
                    /* Level2 S-Function Block: '<S826>/B_82_302' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 1514, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_h[0] = _rtB->B_97_1511_0[0];
    _rtDW->UnitDelay1_DSTATE_k[0] = _rtB->B_97_1521_0[0];
    _rtDW->UnitDelay_DSTATE_h[1] = _rtB->B_97_1511_0[1];
    _rtDW->UnitDelay1_DSTATE_k[1] = _rtB->B_97_1521_0[1];
    _rtDW->UnitDelay_DSTATE_h[2] = _rtB->B_97_1511_0[2];
    _rtDW->UnitDelay1_DSTATE_k[2] = _rtB->B_97_1521_0[2];
    HoldSine = _rtDW->lastSin_f;
    _rtDW->lastSin_f = _rtDW->lastSin_f * _rtP->P_1652 + _rtDW->lastCos_g * _rtP->P_1651;
    _rtDW->lastCos_g = _rtDW->lastCos_g * _rtP->P_1652 - HoldSine * _rtP->P_1651;
    _rtDW->Integ4_SYSTEM_ENABLE_lp = 0U;
    _rtDW->Integ4_DSTATE_hb[0] = _rtP->P_1655 * _rtB->B_97_1523_0[0] + _rtB->B_97_1524_0[0];
    _rtDW->Integ4_DSTATE_hb[1] = _rtP->P_1655 * _rtB->B_97_1523_0[1] + _rtB->B_97_1524_0[1];
    _rtDW->Integ4_DSTATE_hb[2] = _rtP->P_1655 * _rtB->B_97_1523_0[2] + _rtB->B_97_1524_0[2];
                    /* Level2 S-Function Block: '<S824>/B_82_304' (sfun_discreteVariableDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 1526, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_l3k[0] = _rtB->B_97_1523_0[0];
    _rtDW->UnitDelay1_DSTATE_d[0] = _rtB->B_97_1533_0[0];
    _rtDW->UnitDelay_DSTATE_l3k[1] = _rtB->B_97_1523_0[1];
    _rtDW->UnitDelay1_DSTATE_d[1] = _rtB->B_97_1533_0[1];
    _rtDW->UnitDelay_DSTATE_l3k[2] = _rtB->B_97_1523_0[2];
    _rtDW->UnitDelay1_DSTATE_d[2] = _rtB->B_97_1533_0[2];
}
            {
    real_T **uBuffer = (real_T**)&_rtDW->T_PWORK_dy.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T_PWORK_dy.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T_IWORK_i.Head = ((_rtDW->T_IWORK_i.Head < (_rtDW->T_IWORK_i.CircularBufSize-1)) ? (_rtDW->T_IWORK_i.Head+1) : 0);

      if (_rtDW->T_IWORK_i.Head == _rtDW->T_IWORK_i.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T_IWORK_i.CircularBufSize, &_rtDW->T_IWORK_i.Tail, &_rtDW->T_IWORK_i.Head, &_rtDW->T_IWORK_i.Last, simTime - _rtP->P_1674, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T_IWORK_i.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T_IWORK_i.Head] = simTime;
        (*uBuffer)[_rtDW->T_IWORK_i.Head] = _rtB->B_97_1545_0;    
  }
             {
    real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK_ck.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK_ck.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T1_IWORK_g.Head = ((_rtDW->T1_IWORK_g.Head < (_rtDW->T1_IWORK_g.CircularBufSize-1)) ? (_rtDW->T1_IWORK_g.Head+1) : 0);

      if (_rtDW->T1_IWORK_g.Head == _rtDW->T1_IWORK_g.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T1_IWORK_g.CircularBufSize, &_rtDW->T1_IWORK_g.Tail, &_rtDW->T1_IWORK_g.Head, &_rtDW->T1_IWORK_g.Last, simTime - _rtP->P_1677, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T1_IWORK_g.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T1_IWORK_g.Head] = simTime;
        (*uBuffer)[_rtDW->T1_IWORK_g.Head] = _rtB->B_97_1548_0;    
  }
             {
    real_T **uBuffer = (real_T**)&_rtDW->T_PWORK_p.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T_PWORK_p.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T_IWORK_i5.Head = ((_rtDW->T_IWORK_i5.Head < (_rtDW->T_IWORK_i5.CircularBufSize-1)) ? (_rtDW->T_IWORK_i5.Head+1) : 0);

      if (_rtDW->T_IWORK_i5.Head == _rtDW->T_IWORK_i5.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T_IWORK_i5.CircularBufSize, &_rtDW->T_IWORK_i5.Tail, &_rtDW->T_IWORK_i5.Head, &_rtDW->T_IWORK_i5.Last, simTime - _rtP->P_1680, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T_IWORK_i5.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T_IWORK_i5.Head] = simTime;
        (*uBuffer)[_rtDW->T_IWORK_i5.Head] = _rtB->B_97_1552_0;    
  }
             {
    real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK_ph.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK_ph.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T1_IWORK_i.Head = ((_rtDW->T1_IWORK_i.Head < (_rtDW->T1_IWORK_i.CircularBufSize-1)) ? (_rtDW->T1_IWORK_i.Head+1) : 0);

      if (_rtDW->T1_IWORK_i.Head == _rtDW->T1_IWORK_i.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T1_IWORK_i.CircularBufSize, &_rtDW->T1_IWORK_i.Tail, &_rtDW->T1_IWORK_i.Head, &_rtDW->T1_IWORK_i.Last, simTime - _rtP->P_1683, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T1_IWORK_i.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T1_IWORK_i.Head] = simTime;
        (*uBuffer)[_rtDW->T1_IWORK_i.Head] = _rtB->B_97_1555_0;    
  }
             {
    real_T **uBuffer = (real_T**)&_rtDW->T_PWORK_bl.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T_PWORK_bl.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T_IWORK_by.Head = ((_rtDW->T_IWORK_by.Head < (_rtDW->T_IWORK_by.CircularBufSize-1)) ? (_rtDW->T_IWORK_by.Head+1) : 0);

      if (_rtDW->T_IWORK_by.Head == _rtDW->T_IWORK_by.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T_IWORK_by.CircularBufSize, &_rtDW->T_IWORK_by.Tail, &_rtDW->T_IWORK_by.Head, &_rtDW->T_IWORK_by.Last, simTime - _rtP->P_1686, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T_IWORK_by.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T_IWORK_by.Head] = simTime;
        (*uBuffer)[_rtDW->T_IWORK_by.Head] = _rtB->B_97_1559_0;    
  }
             {
    real_T **uBuffer = (real_T**)&_rtDW->T1_PWORK_oh.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->T1_PWORK_oh.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->T1_IWORK_p.Head = ((_rtDW->T1_IWORK_p.Head < (_rtDW->T1_IWORK_p.CircularBufSize-1)) ? (_rtDW->T1_IWORK_p.Head+1) : 0);

      if (_rtDW->T1_IWORK_p.Head == _rtDW->T1_IWORK_p.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->T1_IWORK_p.CircularBufSize, &_rtDW->T1_IWORK_p.Tail, &_rtDW->T1_IWORK_p.Head, &_rtDW->T1_IWORK_p.Last, simTime - _rtP->P_1689, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->T1_IWORK_p.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->T1_IWORK_p.Head] = simTime;
        (*uBuffer)[_rtDW->T1_IWORK_p.Head] = _rtB->B_97_1562_0;    
  }
 i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtDW->UnitDelay_DSTATE_jz = _rtB->B_97_1646_0;
    _rtDW->DiscreteTimeIntegrator_DSTATE_om += _rtP->P_1712 * _rtB->B_97_1640_0;
    _rtDW->Integ4_SYSTEM_ENABLE_b = 0U;
    _rtDW->Integ4_DSTATE_m2 = _rtP->P_1717 * _rtB->B_97_1616_0 + _rtB->B_97_1617_0;
                    /* Level2 S-Function Block: '<S812>/B_82_319' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 97, 1622, SS_CALL_MDL_UPDATE);

     _rtDW->UnitDelay_DSTATE_o0 = _rtB->B_97_1616_0;
    if (_rtDW->AutomaticGainControl_MODE_o) {
        i = ssIsSampleHit(S, 2, 0);
        if (i != 0) {
            _rtDW->Integ4_SYSTEM_ENABLE_au = 0U;
            _rtDW->Integ4_DSTATE_fi = _rtP->P_462 * _rtB->B_77_10_0 + _rtB->B_77_11_0;
                    /* Level2 S-Function Block: '<S807>/B_66_0' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 77, 16, SS_CALL_MDL_UPDATE);

             _rtDW->UnitDelay_DSTATE_hk = _rtB->B_77_10_0;
            _rtDW->Integ4_SYSTEM_ENABLE_d = 0U;
            _rtDW->Integ4_DSTATE_ez = _rtP->P_473 * _rtB->B_77_22_0 + _rtB->B_77_23_0;
                    /* Level2 S-Function Block: '<S809>/B_66_2' (discreteVariableTransportDelay) */
          /* Call into Simulink for MEX-version of S-function */
          ssCallAccelRunBlock(S, 77, 28, SS_CALL_MDL_UPDATE);

             _rtDW->UnitDelay_DSTATE_fvk = _rtB->B_77_22_0;
        }
    }
    _rtDW->DiscreteTimeIntegrator_DSTATE_k += _rtP->P_1730 * _rtB->B_97_1648_0;
    if (_rtDW->DiscreteTimeIntegrator_DSTATE_k >= _rtP->P_1732) {
        _rtDW->DiscreteTimeIntegrator_DSTATE_k = _rtP->P_1732;
    } else {
        if (_rtDW->DiscreteTimeIntegrator_DSTATE_k <= _rtP->P_1733) {
            _rtDW->DiscreteTimeIntegrator_DSTATE_k = _rtP->P_1733;
        }
    }
    _rtDW->UD_DSTATE_k = _rtB->B_97_1636_0;
    _rtDW->UnitDelay_DSTATE_d0 = _rtB->B_97_1645_0;
              {
          real_T xnew[2];
        xnew[0] = (_rtP->P_1743[0])*_rtDW->DiscreteStateSpace_DSTATE_e[0] + (_rtP->P_1743[1])*_rtDW->DiscreteStateSpace_DSTATE_e[1];
    xnew[0] += (_rtP->P_1744[0])*_rtB->B_97_1645_0;

    xnew[1] = (_rtP->P_1743[2])*_rtDW->DiscreteStateSpace_DSTATE_e[0] + (_rtP->P_1743[3])*_rtDW->DiscreteStateSpace_DSTATE_e[1];
    xnew[1] += (_rtP->P_1744[1])*_rtB->B_97_1645_0;
    (void) memcpy(&_rtDW->DiscreteStateSpace_DSTATE_e[0], xnew,
sizeof(real_T)*2);
    }

 }
            {
    real_T **uBuffer = (real_T**)&_rtDW->ENGINETd_PWORK_d.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&_rtDW->ENGINETd_PWORK_d.TUbufferPtrs[1];
       
    real_T simTime   = ssGetT(S);
    
    
        

      
      _rtDW->ENGINETd_IWORK_c.Head = ((_rtDW->ENGINETd_IWORK_c.Head < (_rtDW->ENGINETd_IWORK_c.CircularBufSize-1)) ? (_rtDW->ENGINETd_IWORK_c.Head+1) : 0);

      if (_rtDW->ENGINETd_IWORK_c.Head == _rtDW->ENGINETd_IWORK_c.Tail) {
        
          
          
	  if (!Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_rt_TDelayUpdateTailOrGrowBuf( &_rtDW->ENGINETd_IWORK_c.CircularBufSize, &_rtDW->ENGINETd_IWORK_c.Tail, &_rtDW->ENGINETd_IWORK_c.Head, &_rtDW->ENGINETd_IWORK_c.Last, simTime - _rtP->P_1750, tBuffer, uBuffer, (NULL), (boolean_T)0, false, &_rtDW->ENGINETd_IWORK_c.MaxNewBufSize)) {
            ssSetErrorStatus(S, "tdelay memory allocation error");
            return;
          }
      }
      
      (*tBuffer)[_rtDW->ENGINETd_IWORK_c.Head] = simTime;
        (*uBuffer)[_rtDW->ENGINETd_IWORK_c.Head] = _rtB->B_97_1796_0;    
  }
 i = ssIsSampleHit(S, 2, 0);
if (i != 0) {
    _rtDW->DiscreteTimeIntegrator2_DSTATE_j += _rtP->P_1755 * _rtB->B_97_1811_0;
    _rtDW->UnitDelay2_DSTATE_h = _rtB->B_97_1791_0;
              {
          real_T xnew[1];
        xnew[0] = _rtP->P_1776*_rtDW->DiscreteStateSpace_DSTATE_ax;
    xnew[0] += _rtP->P_1777*_rtB->B_97_1763_0;
    (void) memcpy(&_rtDW->DiscreteStateSpace_DSTATE_ax, xnew,
sizeof(real_T)*1);
    }

     _rtDW->UnitDelay1_DSTATE_fx = _rtB->B_97_1793_0;
              {
          real_T xnew[1];
        xnew[0] = _rtP->P_1788*_rtDW->DiscreteStateSpace_DSTATE_m5;
    xnew[0] += _rtP->P_1789*_rtB->B_97_1768_0;
    (void) memcpy(&_rtDW->DiscreteStateSpace_DSTATE_m5, xnew,
sizeof(real_T)*1);
    }

               {
          real_T xnew[1];
        xnew[0] = (_rtP->P_1799)*_rtDW->DiscreteStateSpace_DSTATE_l;
    xnew[0] += _rtP->P_1800*_rtB->B_97_1790_0;
    (void) memcpy(&_rtDW->DiscreteStateSpace_DSTATE_l, xnew,
sizeof(real_T)*1);
    }

               {
          real_T xnew[1];
        xnew[0] = _rtP->P_1805*_rtDW->DiscreteStateSpace_DSTATE_hk;
    xnew[0] += _rtP->P_1806*_rtB->B_97_1792_0;
    (void) memcpy(&_rtDW->DiscreteStateSpace_DSTATE_hk, xnew,
sizeof(real_T)*1);
    }

 }
/* 0: INTG_NORMAL     1: INTG_LEAVING_UPPER_SAT  2: INTG_LEAVING_LOWER_SAT */
/* 3: INTG_UPPER_SAT  4: INTG_LOWER_SAT */
if (_rtX->Integrator_CSTATE_k == _rtP->P_1814) {
    switch (_rtDW->Integrator_MODE_h) {
      case 3:
        if (_rtB->B_97_1800_0 < 0.0) {
            ssSetBlockStateForSolverChangedAtMajorStep(S);
            _rtDW->Integrator_MODE_h = 1;
        }
        break;
      case 1:
        if (_rtB->B_97_1800_0 >= 0.0) {
            _rtDW->Integrator_MODE_h = 3;
            ssSetBlockStateForSolverChangedAtMajorStep(S);
        }
        break;
      default:
        ssSetBlockStateForSolverChangedAtMajorStep(S);
        if (_rtB->B_97_1800_0 < 0.0) {
            _rtDW->Integrator_MODE_h = 1;
        } else {
            _rtDW->Integrator_MODE_h = 3;
        }
        break;
    }
} else if (_rtX->Integrator_CSTATE_k == _rtP->P_1815) {
    switch (_rtDW->Integrator_MODE_h) {
      case 4:
        if (_rtB->B_97_1800_0 > 0.0) {
            ssSetBlockStateForSolverChangedAtMajorStep(S);
            _rtDW->Integrator_MODE_h = 2;
        }
        break;
      case 2:
        if (_rtB->B_97_1800_0 <= 0.0) {
            _rtDW->Integrator_MODE_h = 4;
            ssSetBlockStateForSolverChangedAtMajorStep(S);
        }
        break;
      default:
        ssSetBlockStateForSolverChangedAtMajorStep(S);
        if (_rtB->B_97_1800_0 > 0.0) {
            _rtDW->Integrator_MODE_h = 2;
        } else {
            _rtDW->Integrator_MODE_h = 4;
        }
        break;
    }
} else {
    _rtDW->Integrator_MODE_h = 0;
}
if (_rtDW->TimeStampA == (rtInf)) {
    _rtDW->TimeStampA = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA;
} else if (_rtDW->TimeStampB == (rtInf)) {
    _rtDW->TimeStampB = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB;
} else if (_rtDW->TimeStampA < _rtDW->TimeStampB) {
    _rtDW->TimeStampA = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA;
} else {
    _rtDW->TimeStampB = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB;
}
*lastU = _rtB->B_97_1905_0;
if (_rtDW->TimeStampA_c == (rtInf)) {
    _rtDW->TimeStampA_c = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA_d;
} else if (_rtDW->TimeStampB_a == (rtInf)) {
    _rtDW->TimeStampB_a = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB_d;
} else if (_rtDW->TimeStampA_c < _rtDW->TimeStampB_a) {
    _rtDW->TimeStampA_c = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA_d;
} else {
    _rtDW->TimeStampB_a = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB_d;
}
*lastU = _rtB->B_97_1917_0;
if (_rtDW->TimeStampA_cb == (rtInf)) {
    _rtDW->TimeStampA_cb = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA_d2;
} else if (_rtDW->TimeStampB_ax == (rtInf)) {
    _rtDW->TimeStampB_ax = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB_m;
} else if (_rtDW->TimeStampA_cb < _rtDW->TimeStampB_ax) {
    _rtDW->TimeStampA_cb = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA_d2;
} else {
    _rtDW->TimeStampB_ax = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB_m;
}
*lastU = _rtB->B_97_1929_0;
if (_rtDW->TimeStampA_co == (rtInf)) {
    _rtDW->TimeStampA_co = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA_l;
} else if (_rtDW->TimeStampB_f == (rtInf)) {
    _rtDW->TimeStampB_f = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB_i;
} else if (_rtDW->TimeStampA_co < _rtDW->TimeStampB_f) {
    _rtDW->TimeStampA_co = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA_l;
} else {
    _rtDW->TimeStampB_f = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB_i;
}
*lastU = _rtB->B_97_1954_0;
if (_rtDW->TimeStampA_k == (rtInf)) {
    _rtDW->TimeStampA_k = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA_m;
} else if (_rtDW->TimeStampB_g == (rtInf)) {
    _rtDW->TimeStampB_g = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB_n;
} else if (_rtDW->TimeStampA_k < _rtDW->TimeStampB_g) {
    _rtDW->TimeStampA_k = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA_m;
} else {
    _rtDW->TimeStampB_g = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB_n;
}
*lastU = _rtB->B_97_1966_0;
if (_rtDW->TimeStampA_n == (rtInf)) {
    _rtDW->TimeStampA_n = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA_n;
} else if (_rtDW->TimeStampB_e == (rtInf)) {
    _rtDW->TimeStampB_e = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB_nm;
} else if (_rtDW->TimeStampA_n < _rtDW->TimeStampB_e) {
    _rtDW->TimeStampA_n = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA_n;
} else {
    _rtDW->TimeStampB_e = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB_nm;
}
*lastU = _rtB->B_97_1978_0;
if (_rtDW->TimeStampA_k3 == (rtInf)) {
    _rtDW->TimeStampA_k3 = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA_mm;
} else if (_rtDW->TimeStampB_j == (rtInf)) {
    _rtDW->TimeStampB_j = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB_l;
} else if (_rtDW->TimeStampA_k3 < _rtDW->TimeStampB_j) {
    _rtDW->TimeStampA_k3 = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA_mm;
} else {
    _rtDW->TimeStampB_j = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB_l;
}
*lastU = _rtB->B_97_2003_0;
if (_rtDW->TimeStampA_l == (rtInf)) {
    _rtDW->TimeStampA_l = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA_h;
} else if (_rtDW->TimeStampB_l == (rtInf)) {
    _rtDW->TimeStampB_l = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB_b;
} else if (_rtDW->TimeStampA_l < _rtDW->TimeStampB_l) {
    _rtDW->TimeStampA_l = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA_h;
} else {
    _rtDW->TimeStampB_l = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB_b;
}
*lastU = _rtB->B_97_2015_0;
if (_rtDW->TimeStampA_cm == (rtInf)) {
    _rtDW->TimeStampA_cm = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA_n3;
} else if (_rtDW->TimeStampB_gf == (rtInf)) {
    _rtDW->TimeStampB_gf = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB_g;
} else if (_rtDW->TimeStampA_cm < _rtDW->TimeStampB_gf) {
    _rtDW->TimeStampA_cm = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA_n3;
} else {
    _rtDW->TimeStampB_gf = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB_g;
}
*lastU = _rtB->B_97_2027_0;
if (_rtDW->TimeStampA_f == (rtInf)) {
    _rtDW->TimeStampA_f = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA_dj;
} else if (_rtDW->TimeStampB_gp == (rtInf)) {
    _rtDW->TimeStampB_gp = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB_o;
} else if (_rtDW->TimeStampA_f < _rtDW->TimeStampB_gp) {
    _rtDW->TimeStampA_f = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA_dj;
} else {
    _rtDW->TimeStampB_gp = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB_o;
}
*lastU = _rtB->B_97_2052_0;
if (_rtDW->TimeStampA_d == (rtInf)) {
    _rtDW->TimeStampA_d = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA_dy;
} else if (_rtDW->TimeStampB_eg == (rtInf)) {
    _rtDW->TimeStampB_eg = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB_a;
} else if (_rtDW->TimeStampA_d < _rtDW->TimeStampB_eg) {
    _rtDW->TimeStampA_d = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA_dy;
} else {
    _rtDW->TimeStampB_eg = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB_a;
}
*lastU = _rtB->B_97_2064_0;
if (_rtDW->TimeStampA_b == (rtInf)) {
    _rtDW->TimeStampA_b = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA_f;
} else if (_rtDW->TimeStampB_k == (rtInf)) {
    _rtDW->TimeStampB_k = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB_h;
} else if (_rtDW->TimeStampA_b < _rtDW->TimeStampB_k) {
    _rtDW->TimeStampA_b = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeA_f;
} else {
    _rtDW->TimeStampB_k = ssGetT(S);
    lastU = &_rtDW->LastUAtTimeB_h;
}
*lastU = _rtB->B_97_2076_0;







    
  

                  UNUSED_PARAMETER(tid);



  

            }
      
  






    
  

          /* Update for root system: '<Root>' */
      
            #define MDL_UPDATE
    
  
  
      
  
     static  void mdlUpdateTID5(SimStruct *S, int_T tid)
  {
  
      
        
      
    
  

            
      
    
  

                
  

              UNUSED_PARAMETER(tid);



  

            }
      
  








	      /* Derivatives for root system: '<Root>' */
              #define MDL_DERIVATIVES
    
  
  
      
  
     static  void mdlDerivatives(SimStruct *S)
  {
  
  
      
  
      
  
B_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *_rtB;
P_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *_rtP;
X_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *_rtX;
XDot_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *_rtXdot;
XDis_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *_rtXdis;
DW_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *_rtDW;

    
  

              
    
    
  

                    
                _rtDW = ((DW_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetRootDWork(S));
_rtXdis = ((XDis_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetContStateDisabled(S));
_rtXdot = ((XDot_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetdX(S));
_rtX = ((X_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetContStates(S));
_rtP = ((P_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetModelRtp(S));
_rtB = ((B_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) _ssGetModelBlockIO(S));
_rtXdot->integ1_CSTATE = _rtB->B_97_430_0;
_rtXdot->Integ2_CSTATE = _rtB->B_97_432_0;
_rtXdot->integ1_CSTATE_n = _rtB->B_97_434_0;
_rtXdot->Integ2_CSTATE_p = _rtB->B_97_435_0;
_rtXdot->integ1_CSTATE_o = _rtB->B_97_437_0;
_rtXdot->Integ2_CSTATE_h = _rtB->B_97_438_0;
_rtXdot->TransferFcn1_CSTATE = 0.0;
_rtXdot->TransferFcn1_CSTATE += _rtP->P_809 * _rtX->TransferFcn1_CSTATE;
_rtXdot->TransferFcn1_CSTATE += _rtB->B_97_458_0;
if (_rtDW->GridSupportingasCurrentSourceGridFeeding_MODE) {
    _rtXdot->TransferFcn1_CSTATE_p = 0.0;
    _rtXdot->TransferFcn1_CSTATE_p += _rtP->P_58 * _rtX->TransferFcn1_CSTATE_p;
    _rtXdot->TransferFcn1_CSTATE_p += _rtB->B_97_363_0;
    _rtXdot->Integrator_CSTATE_e = _rtB->B_27_32_0;
    _rtXdot->Filter_CSTATE_pr = _rtB->B_27_10_0;
    _rtXdot->TransferFcn2_CSTATE_f = 0.0;
    _rtXdot->TransferFcn2_CSTATE_f += _rtP->P_67 * _rtX->TransferFcn2_CSTATE_f;
    _rtXdot->TransferFcn2_CSTATE_f += _rtB->B_27_37_0;
    _rtXdot->Integrator_CSTATE_o = _rtB->B_27_31_0;
    _rtXdot->Filter_CSTATE_d = _rtB->B_27_17_0;
    _rtXdot->TransferFcn4_CSTATE_p = 0.0;
    _rtXdot->TransferFcn4_CSTATE_p += _rtP->P_74 * _rtX->TransferFcn4_CSTATE_p;
    _rtXdot->TransferFcn4_CSTATE_p += _rtB->B_27_35_0;
    _rtXdot->Integrator_CSTATE_p = _rtB->B_27_33_0;
    _rtXdot->Filter_CSTATE_j = _rtB->B_27_24_0;
    _rtXdot->TransferFcn4_CSTATE_d = 0.0;
    _rtXdot->TransferFcn4_CSTATE_d += _rtP->P_89 * _rtX->TransferFcn4_CSTATE_d;
    _rtXdot->TransferFcn4_CSTATE_d += _rtB->B_27_241_0;
    _rtXdot->Integrator_CSTATE_pv = _rtB->B_27_239_0;
    _rtXdot->Filter_CSTATE_l = _rtB->B_27_45_0;
    _rtXdot->Integrator_CSTATE_ou = _rtB->B_27_238_0;
    _rtXdot->Filter_CSTATE_f = _rtB->B_27_55_0;
    _rtXdot->TransferFcn2_CSTATE_b = 0.0;
    _rtXdot->TransferFcn2_CSTATE_b += _rtP->P_102 * _rtX->TransferFcn2_CSTATE_b;
    _rtXdot->TransferFcn2_CSTATE_b += _rtB->B_27_243_0;
    _rtXdot->Integrator_CSTATE_cm = _rtB->B_27_237_0;
    _rtXdot->Filter_CSTATE_n = _rtB->B_27_62_0;
} else {
      {
        real_T *dx;
        int_T   i;
          
  

  
      dx = &(((XDot_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetdX(S))->TransferFcn1_CSTATE_p);
    for (i=0; i < 17; i++) {
    dx[i] = 0.0;
    }
  
      }
}
if (_rtDW->Gridforming_MODE) {
    _rtXdot->TransferFcn5_CSTATE = 0.0;
    _rtXdot->TransferFcn5_CSTATE += _rtP->P_267 * _rtX->TransferFcn5_CSTATE;
    _rtXdot->TransferFcn5_CSTATE += _rtB->B_36_102_0;
    _rtXdot->TransferFcn4_CSTATE_f = 0.0;
    _rtXdot->TransferFcn4_CSTATE_f += _rtP->P_279 * _rtX->TransferFcn4_CSTATE_f;
    _rtXdot->TransferFcn4_CSTATE_f += _rtB->B_36_115_0[0];
    _rtXdot->TransferFcn6_CSTATE = 0.0;
    _rtXdot->TransferFcn6_CSTATE += _rtP->P_281 * _rtX->TransferFcn6_CSTATE;
    _rtXdot->TransferFcn6_CSTATE += _rtB->B_36_115_0[1];
    _rtXdot->TransferFcn7_CSTATE = 0.0;
    _rtXdot->TransferFcn7_CSTATE += _rtP->P_286 * _rtX->TransferFcn7_CSTATE;
    _rtXdot->TransferFcn7_CSTATE += _rtB->B_36_115_0[2];
} else {
      {
        real_T *dx;
        int_T   i;
          
  

  
      dx = &(((XDot_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetdX(S))->TransferFcn5_CSTATE);
    for (i=0; i < 4; i++) {
    dx[i] = 0.0;
    }
  
      }
}
_rtXdot->integ1_CSTATE_nd = _rtB->B_97_816_0;
_rtXdot->Integ2_CSTATE_i = _rtB->B_97_818_0;
_rtXdot->integ1_CSTATE_d = _rtB->B_97_820_0;
_rtXdot->Integ2_CSTATE_pa = _rtB->B_97_821_0;
_rtXdot->integ1_CSTATE_e = _rtB->B_97_823_0;
_rtXdot->Integ2_CSTATE_e = _rtB->B_97_824_0;
_rtXdot->TransferFcn1_CSTATE_k = 0.0;
_rtXdot->TransferFcn1_CSTATE_k += _rtP->P_1080 * _rtX->TransferFcn1_CSTATE_k;
_rtXdot->TransferFcn1_CSTATE_k += _rtB->B_97_844_0;
_rtXdot->TransferFcn4_CSTATE = 0.0;
_rtXdot->TransferFcn4_CSTATE += _rtP->P_1159 * _rtX->TransferFcn4_CSTATE;
_rtXdot->TransferFcn4_CSTATE += _rtB->B_97_953_0;
_rtXdot->Integrator_CSTATE = _rtB->B_97_951_0;
_rtXdot->Filter_CSTATE = _rtB->B_97_888_0;
_rtXdot->TransferFcn1_CSTATE_m = 0.0;
_rtXdot->TransferFcn1_CSTATE_m += _rtP->P_1165 * _rtX->TransferFcn1_CSTATE_m;
_rtXdot->TransferFcn1_CSTATE_m += _rtB->B_97_761_0;
_rtXdot->Integrator_CSTATE_c = _rtB->B_97_950_0;
_rtXdot->Filter_CSTATE_e = _rtB->B_97_899_0;
_rtXdot->TransferFcn2_CSTATE = 0.0;
_rtXdot->TransferFcn2_CSTATE += _rtP->P_1174 * _rtX->TransferFcn2_CSTATE;
_rtXdot->TransferFcn2_CSTATE += _rtB->B_97_955_0;
_rtXdot->Integrator_CSTATE_a = _rtB->B_97_949_0;
_rtXdot->Filter_CSTATE_g = _rtB->B_97_906_0;
_rtXdot->TransferFcn4_CSTATE_k = 0.0;
_rtXdot->TransferFcn4_CSTATE_k += _rtP->P_1184 * _rtX->TransferFcn4_CSTATE_k;
_rtXdot->TransferFcn4_CSTATE_k += _rtB->B_97_1126_0;
_rtXdot->Integrator_CSTATE_h = _rtB->B_97_1124_0;
_rtXdot->Filter_CSTATE_c = _rtB->B_97_920_0;
_rtXdot->Integrator_CSTATE_g = _rtB->B_97_1123_0;
_rtXdot->Filter_CSTATE_b = _rtB->B_97_930_0;
_rtXdot->TransferFcn2_CSTATE_d = 0.0;
_rtXdot->TransferFcn2_CSTATE_d += _rtP->P_1197 * _rtX->TransferFcn2_CSTATE_d;
_rtXdot->TransferFcn2_CSTATE_d += _rtB->B_97_1128_0;
_rtXdot->Integrator_CSTATE_h0 = _rtB->B_97_1122_0;
_rtXdot->Filter_CSTATE_p = _rtB->B_97_937_0;
_rtXdot->CONTROLSYSTEM_CSTATE[0] = 0.0;
_rtXdot->CONTROLSYSTEM_CSTATE[0] += _rtP->P_1374[0] * _rtX->CONTROLSYSTEM_CSTATE[0];
_rtXdot->CONTROLSYSTEM_CSTATE[1] = 0.0;
_rtXdot->CONTROLSYSTEM_CSTATE[0] += _rtP->P_1374[1] * _rtX->CONTROLSYSTEM_CSTATE[1];
_rtXdot->CONTROLSYSTEM_CSTATE[1] += _rtX->CONTROLSYSTEM_CSTATE[0];
_rtXdot->CONTROLSYSTEM_CSTATE[0] += _rtB->B_97_1246_0;
/* 0: INTG_NORMAL     1: INTG_LEAVING_UPPER_SAT  2: INTG_LEAVING_LOWER_SAT */
/* 3: INTG_UPPER_SAT  4: INTG_LOWER_SAT */
if ((_rtDW->Integrator_MODE != 3) && (_rtDW->Integrator_MODE != 4)) {
    _rtXdot->Integrator_CSTATE_g1 = _rtB->B_97_1248_0;
    _rtXdis->Integrator_CSTATE_g1 = false;
} else {
    /* in saturation */
    _rtXdot->Integrator_CSTATE_g1 = 0.0;
    if ((_rtDW->Integrator_MODE == 3) || (_rtDW->Integrator_MODE == 4)) {
        _rtXdis->Integrator_CSTATE_g1 = true;
    }
}
_rtXdot->TF1_CSTATE = 0.0;
_rtXdot->TF1_CSTATE += _rtP->P_1390 * _rtX->TF1_CSTATE;
_rtXdot->TF1_CSTATE += _rtB->B_97_1231_0;
_rtXdot->TF2_CSTATE = 0.0;
_rtXdot->TF2_CSTATE += _rtP->P_1393 * _rtX->TF2_CSTATE;
_rtXdot->TF2_CSTATE += _rtB->B_97_1247_0;
_rtXdot->integ1_CSTATE_oz = _rtB->B_97_1379_0;
_rtXdot->Integ2_CSTATE_k = _rtB->B_97_1381_0;
_rtXdot->integ1_CSTATE_nb = _rtB->B_97_1383_0;
_rtXdot->Integ2_CSTATE_ik = _rtB->B_97_1384_0;
_rtXdot->integ1_CSTATE_g = _rtB->B_97_1386_0;
_rtXdot->Integ2_CSTATE_o = _rtB->B_97_1387_0;
_rtXdot->TransferFcn1_CSTATE_f = 0.0;
_rtXdot->TransferFcn1_CSTATE_f += _rtP->P_1520 * _rtX->TransferFcn1_CSTATE_f;
_rtXdot->TransferFcn1_CSTATE_f += _rtB->B_97_1407_0;
_rtXdot->integ1_CSTATE_do = _rtB->B_97_1590_0;
_rtXdot->Integ2_CSTATE_oo = _rtB->B_97_1592_0;
_rtXdot->integ1_CSTATE_m = _rtB->B_97_1594_0;
_rtXdot->Integ2_CSTATE_ot = _rtB->B_97_1595_0;
_rtXdot->integ1_CSTATE_i = _rtB->B_97_1597_0;
_rtXdot->Integ2_CSTATE_j = _rtB->B_97_1598_0;
_rtXdot->CONTROLSYSTEM_CSTATE_b[0] = 0.0;
_rtXdot->CONTROLSYSTEM_CSTATE_b[0] += _rtP->P_1810[0] * _rtX->CONTROLSYSTEM_CSTATE_b[0];
_rtXdot->CONTROLSYSTEM_CSTATE_b[1] = 0.0;
_rtXdot->CONTROLSYSTEM_CSTATE_b[0] += _rtP->P_1810[1] * _rtX->CONTROLSYSTEM_CSTATE_b[1];
_rtXdot->CONTROLSYSTEM_CSTATE_b[1] += _rtX->CONTROLSYSTEM_CSTATE_b[0];
_rtXdot->CONTROLSYSTEM_CSTATE_b[0] += _rtB->B_97_1798_0;
/* 0: INTG_NORMAL     1: INTG_LEAVING_UPPER_SAT  2: INTG_LEAVING_LOWER_SAT */
/* 3: INTG_UPPER_SAT  4: INTG_LOWER_SAT */
if ((_rtDW->Integrator_MODE_h != 3) && (_rtDW->Integrator_MODE_h != 4)) {
    _rtXdot->Integrator_CSTATE_k = _rtB->B_97_1800_0;
    _rtXdis->Integrator_CSTATE_k = false;
} else {
    /* in saturation */
    _rtXdot->Integrator_CSTATE_k = 0.0;
    if ((_rtDW->Integrator_MODE_h == 3) || (_rtDW->Integrator_MODE_h == 4)) {
        _rtXdis->Integrator_CSTATE_k = true;
    }
}
_rtXdot->TF1_CSTATE_p = 0.0;
_rtXdot->TF1_CSTATE_p += _rtP->P_1817 * _rtX->TF1_CSTATE_p;
_rtXdot->TF1_CSTATE_p += _rtB->B_97_1795_0;
_rtXdot->TF2_CSTATE_m = 0.0;
_rtXdot->TF2_CSTATE_m += _rtP->P_1820 * _rtX->TF2_CSTATE_m;
_rtXdot->TF2_CSTATE_m += _rtB->B_97_1799_0;






    
  

      
      
      }
      
  



	
	      /* ZeroCrossings for root system: '<Root>' */
            #define MDL_ZERO_CROSSINGS
    
  
  
      
  
     static  void mdlZeroCrossings(SimStruct *S)
  {
  
  
  
        
  
boolean_T anyStateSaturated;
B_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *_rtB;
P_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *_rtP;
X_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *_rtX;
ZCV_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *_rtZCSV;
DW_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *_rtDW;

    
  

                
  

                    
                _rtDW = ((DW_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetRootDWork(S));
_rtZCSV = ((ZCV_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetSolverZcSignalVector(S));
_rtX = ((X_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetContStates(S));
_rtP = ((P_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetModelRtp(S));
_rtB = ((B_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) _ssGetModelBlockIO(S));
_rtZCSV->GridSupportingasCurrentSourceGridFeeding_Enable_ZC = _rtB->B_97_407_0;
if (_rtDW->GridSupportingasCurrentSourceGridFeeding_MODE) {
    _rtZCSV->Saturation_UprLim_ZC_p = _rtB->B_27_72_0 - _rtP->P_113;
    _rtZCSV->Saturation_LwrLim_ZC_e = _rtB->B_27_72_0 - _rtP->P_114;
    _rtZCSV->Saturation1_UprLim_ZC_e = _rtB->B_27_74_0 - _rtP->P_115;
    _rtZCSV->Saturation1_LwrLim_ZC_l = _rtB->B_27_74_0 - _rtP->P_116;
} else {
    
    {
        
  
  
      real_T* zcsv = &(((ZCV_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetSolverZcSignalVector(S))->Saturation_UprLim_ZC_p);
    int_T i;
    for (i=0; i < 4; i++) {
      zcsv[i] = 0.0;
    }

    }
    
}
_rtZCSV->If1_IfInput_ZC = 0.0;
if (_rtB->B_97_405_0 > 0.5) {
    _rtZCSV->If1_IfInput_ZC = 1.0;
}
_rtZCSV->RelationalOperator_RelopInput_ZC[0] = _rtB->B_97_361_0[1] - _rtB->B_97_425_0;
_rtZCSV->RelationalOperator_RelopInput_ZC[1] = _rtB->B_97_361_0[2] - _rtB->B_97_425_0;
_rtZCSV->RelationalOperator_RelopInput_ZC_o[0] = _rtB->B_97_759_0[1] - _rtB->B_97_811_0;
_rtZCSV->RelationalOperator_RelopInput_ZC_o[1] = _rtB->B_97_759_0[2] - _rtB->B_97_811_0;
_rtZCSV->Saturation_UprLim_ZC = _rtB->B_97_957_0 - _rtP->P_1215;
_rtZCSV->Saturation_LwrLim_ZC = _rtB->B_97_957_0 - _rtP->P_1216;
_rtZCSV->Saturation1_UprLim_ZC = _rtB->B_97_960_0 - _rtP->P_1218;
_rtZCSV->Saturation1_LwrLim_ZC = _rtB->B_97_960_0 - _rtP->P_1219;
/* zero crossings for entering into limited region */
/* 0: INTG_NORMAL     1: INTG_LEAVING_UPPER_SAT  2: INTG_LEAVING_LOWER_SAT */
/* 3: INTG_UPPER_SAT  4: INTG_LOWER_SAT */
if ((_rtDW->Integrator_MODE == 1) && (_rtX->Integrator_CSTATE_g1 >= _rtP->P_1380)) {
    _rtZCSV->Integrator_IntgUpLimit_ZC = 0.0;
} else {
    _rtZCSV->Integrator_IntgUpLimit_ZC = _rtX->Integrator_CSTATE_g1 - _rtP->P_1380;
}
if ((_rtDW->Integrator_MODE == 2) && (_rtX->Integrator_CSTATE_g1 <= _rtP->P_1381)) {
    _rtZCSV->Integrator_IntgLoLimit_ZC = 0.0;
} else {
    _rtZCSV->Integrator_IntgLoLimit_ZC = _rtX->Integrator_CSTATE_g1 - _rtP->P_1381;
}
/* zero crossings for leaving limited region */
anyStateSaturated = false;
if ((_rtDW->Integrator_MODE == 3) || (_rtDW->Integrator_MODE == 4)) {
    anyStateSaturated = true;
}
if (anyStateSaturated) {
    _rtZCSV->Integrator_LeaveSaturate_ZC = _rtB->B_97_1248_0;
} else {
    _rtZCSV->Integrator_LeaveSaturate_ZC = 0.0;
}
_rtZCSV->RelationalOperator_RelopInput_ZC_i[0] = _rtB->B_97_1364_0[1] - _rtB->B_97_1374_0;
_rtZCSV->RelationalOperator_RelopInput_ZC_i[1] = _rtB->B_97_1364_0[2] - _rtB->B_97_1374_0;
_rtZCSV->RelationalOperator_RelopInput_ZC_b[0] = _rtB->B_97_1575_0[1] - _rtB->B_97_1585_0;
_rtZCSV->RelationalOperator_RelopInput_ZC_b[1] = _rtB->B_97_1575_0[2] - _rtB->B_97_1585_0;
/* zero crossings for entering into limited region */
/* 0: INTG_NORMAL     1: INTG_LEAVING_UPPER_SAT  2: INTG_LEAVING_LOWER_SAT */
/* 3: INTG_UPPER_SAT  4: INTG_LOWER_SAT */
if ((_rtDW->Integrator_MODE_h == 1) && (_rtX->Integrator_CSTATE_k >= _rtP->P_1814)) {
    _rtZCSV->Integrator_IntgUpLimit_ZC_j = 0.0;
} else {
    _rtZCSV->Integrator_IntgUpLimit_ZC_j = _rtX->Integrator_CSTATE_k - _rtP->P_1814;
}
if ((_rtDW->Integrator_MODE_h == 2) && (_rtX->Integrator_CSTATE_k <= _rtP->P_1815)) {
    _rtZCSV->Integrator_IntgLoLimit_ZC_o = 0.0;
} else {
    _rtZCSV->Integrator_IntgLoLimit_ZC_o = _rtX->Integrator_CSTATE_k - _rtP->P_1815;
}
/* zero crossings for leaving limited region */
anyStateSaturated = false;
if ((_rtDW->Integrator_MODE_h == 3) || (_rtDW->Integrator_MODE_h == 4)) {
    anyStateSaturated = true;
}
if (anyStateSaturated) {
    _rtZCSV->Integrator_LeaveSaturate_ZC_a = _rtB->B_97_1800_0;
} else {
    _rtZCSV->Integrator_LeaveSaturate_ZC_a = 0.0;
}
_rtZCSV->HitCrossing_HitNoOutput_ZC = _rtB->B_97_1907_0 - _rtP->P_1841;
_rtZCSV->HitCrossing_HitNoOutput_ZC_b = _rtB->B_97_1919_0 - _rtP->P_1848;
_rtZCSV->HitCrossing_HitNoOutput_ZC_f = _rtB->B_97_1931_0 - _rtP->P_1855;
_rtZCSV->HitCrossing_HitNoOutput_ZC_h = _rtB->B_97_1956_0 - _rtP->P_1868;
_rtZCSV->HitCrossing_HitNoOutput_ZC_o = _rtB->B_97_1968_0 - _rtP->P_1875;
_rtZCSV->HitCrossing_HitNoOutput_ZC_i = _rtB->B_97_1980_0 - _rtP->P_1882;
_rtZCSV->HitCrossing_HitNoOutput_ZC_f4 = _rtB->B_97_2005_0 - _rtP->P_1895;
_rtZCSV->HitCrossing_HitNoOutput_ZC_h3 = _rtB->B_97_2017_0 - _rtP->P_1902;
_rtZCSV->HitCrossing_HitNoOutput_ZC_l = _rtB->B_97_2029_0 - _rtP->P_1909;
_rtZCSV->HitCrossing_HitNoOutput_ZC_hl = _rtB->B_97_2054_0 - _rtP->P_1922;
_rtZCSV->HitCrossing_HitNoOutput_ZC_n = _rtB->B_97_2066_0 - _rtP->P_1929;
_rtZCSV->HitCrossing_HitNoOutput_ZC_d = _rtB->B_97_2078_0 - _rtP->P_1936;






    
  

          }
      
  



      
/* Function to initialize sizes */
static void mdlInitializeSizes(SimStruct *S)
{

  /* checksum */
  ssSetChecksumVal(S, 0, 1499612481U);
  ssSetChecksumVal(S, 1, 465531497U);
  ssSetChecksumVal(S, 2, 2305553075U);
  ssSetChecksumVal(S, 3, 1004438245U);
  
  { 
     /* First check to see if this is the sample verion of Simulink used 
      * to generate this file. If not force a rebuild and return. */
      mxArray *slVerStructMat = NULL;
      mxArray *slStrMat = mxCreateString("simulink");
      char slVerChar[10];
      int status = mexCallMATLAB(1,&slVerStructMat,
                                 1,&slStrMat,
                                 "ver");
      if(status == 0) {
          mxArray * slVerMat = mxGetField(slVerStructMat,0,"Version");
          if(slVerMat == NULL) {
            status = 1;
          } else {
            status = mxGetString(slVerMat, slVerChar, 10);
          }
      }      
      mxDestroyArray(slStrMat);
      mxDestroyArray(slVerStructMat);
      
      if((status == 1) || (strcmp(slVerChar,"9.3") != 0)) {
          return;
      } 
  }
  
  /* options */
  ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);

    /* Accelerator check memory map size match for DWork */
      if (ssGetSizeofDWork(S) != sizeof(DW_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T)) {
	ssSetErrorStatus(S,"Unexpected error: Internal DWork sizes do "
	"not match for accelerator mex file.");
      }

    /* Accelerator check memory map size match for BlockIO */
    if (ssGetSizeofGlobalBlockIO(S) != sizeof(B_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T)) {
      ssSetErrorStatus(S,"Unexpected error: Internal BlockIO sizes do "
      "not match for accelerator mex file.");
    }
  
  

  
    /* Accelerator check memory map size match for Parameters */
    {
      int ssSizeofParams;
      ssGetSizeofParams(S,&ssSizeofParams);
      if (ssSizeofParams != sizeof(P_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T)) { 
         static char msg[256];
         sprintf(msg,"Unexpected error: Internal Parameters sizes do "
         "not match for accelerator mex file.");
       }
     }
    
    /* model parameters */
    _ssSetModelRtp(S, (real_T *) &Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_rtDefaultP);
    
    /* constant block I/O */
    _ssSetConstBlockIO(S, &Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_rtInvariant);

      
    
  /* non-finites */
  
  rt_InitInfAndNaN(sizeof(real_T));    
    /* non-finite (run-time) assignments */
    
((P_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetModelRtp(S))->P_4 = rtMinusInf;



((P_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetModelRtp(S))->P_590 = rtMinusInf;



((P_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetModelRtp(S))->P_814 = rtInf;



((P_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *) ssGetModelRtp(S))->P_1085 = rtInf;




}

/* mdlInitializeSampleTimes function (used to set up function-call connections) */
static void mdlInitializeSampleTimes(SimStruct *S) { 
          
        slAccRegPrmChangeFcn(S, mdlOutputsTID5);
}

  /* Empty mdlTerminate function (never called) */
  static void mdlTerminate(SimStruct *S) { }

/* MATLAB MEX Glue */
#include "simulink.c"


  

  

  

  
