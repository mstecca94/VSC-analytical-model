

    /*
  * Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_private.h
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


  #ifndef RTW_HEADER_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_private_h_
  #define RTW_HEADER_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_private_h_
  

    

      #include "rtwtypes.h"
  #include "multiword_types.h"

  
  
    #include "Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc.h"


  

  

  

    
      #if !defined(ss_VALIDATE_MEMORY)
  #define ss_VALIDATE_MEMORY(S, ptr)   if(!(ptr)) {\
  ssSetErrorStatus(S, RT_MEMORY_ALLOCATION_ERROR);\
  }
  #endif
  
  #if !defined(rt_FREE)
  #if !defined(_WIN32)
  #define rt_FREE(ptr)   if((ptr) != (NULL)) {\
  free((ptr));\
  (ptr) = (NULL);\
  }
  #else
  /* Visual and other windows compilers declare free without const */
  #define rt_FREE(ptr)   if((ptr) != (NULL)) {\
  free((void *)(ptr));\
  (ptr) = (NULL);\
  }
  #endif
  #endif
  



  

  

  

  

  

  

  

  
          #ifndef __RTW_UTFREE__
      extern void * utMalloc(size_t);
      extern void   utFree(void *); 
      #endif
  

      
    
        

             extern void rt_invd5x5_snf(const real_T u[25], real_T y[25]);
      

 

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
int_T     *maxNewBufSzPtr);

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
;

            
  
        void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_Saturation_Init(SimStruct *S, DW_Saturation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localDW, P_Saturation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localP);
        
  

  
  
        void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_Saturation_Update(SimStruct *S, real_T rtu_Enable, B_Saturation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localB, DW_Saturation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localDW);
        
  

  
  
        void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_Saturation(SimStruct *S, real_T rtu_Enable, const real_T rtu_phi[5], B_Saturation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localB, DW_Saturation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localDW, P_Saturation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localP);
        
  

  
  
        void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_SaturationTID5(SimStruct *S, B_Saturation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localB, P_Saturation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localP);
        
  

  
  
        void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_IfActionSubsystem2_Enable(SimStruct *S);
        
  

  
  
        void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_IfActionSubsystem2_Disable(SimStruct *S);
        
  

  
  
        void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_IfActionSubsystem2(SimStruct *S, real_T rtu_In1, real_T *rty_Out1);
        
  

  
  
        void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_NegSeqComputation_Disable(SimStruct *S, DW_NegSeqComputation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localDW);
        
  

  
  
        void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_NegSeqComputation(SimStruct *S, real_T rtu_Enable, creal_T rtu_In, creal_T rtu_In_h, creal_T rtu_In_l, B_NegSeqComputation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localB, DW_NegSeqComputation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localDW, P_NegSeqComputation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localP);
        
  

  
  
        void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_ZeroSeqComputation_Disable(SimStruct *S, DW_ZeroSeqComputation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localDW);
        
  

  
  
        void Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_ZeroSeqComputation(SimStruct *S, real_T rtu_Enable, creal_T rtu_In, creal_T rtu_In_f, creal_T rtu_In_o, B_ZeroSeqComputation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localB, DW_ZeroSeqComputation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localDW, P_ZeroSeqComputation_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_T *localP);
        
  



  

  

  

  

  

  

  

    #endif /* RTW_HEADER_Improved_With_PV_Disconnection_at_11s_Assignment_1_19052017_acc_private_h_ */
