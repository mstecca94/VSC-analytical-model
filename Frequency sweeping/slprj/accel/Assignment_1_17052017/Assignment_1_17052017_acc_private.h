#include "__cf_Assignment_1_17052017.h"
#ifndef RTW_HEADER_Assignment_1_17052017_acc_private_h_
#define RTW_HEADER_Assignment_1_17052017_acc_private_h_
#include "rtwtypes.h"
#include "multiword_types.h"
#include "Assignment_1_17052017_acc.h"
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
#define rt_FREE(ptr)   if((ptr) != (NULL)) {\
  free((void *)(ptr));\
  (ptr) = (NULL);\
  }
#endif
#endif
#ifndef __RTW_UTFREE__
extern void * utMalloc ( size_t ) ; extern void utFree ( void * ) ;
#endif
extern void rt_invd5x5_snf ( const real_T u [ 25 ] , real_T y [ 25 ] ) ;
boolean_T Assignment_1_17052017_acc_rt_TDelayUpdateTailOrGrowBuf ( int_T *
bufSzPtr , int_T * tailPtr , int_T * headPtr , int_T * lastPtr , real_T
tMinusDelay , real_T * * tBufPtr , real_T * * uBufPtr , real_T * * xBufPtr ,
boolean_T isfixedbuf , boolean_T istransportdelay , int_T * maxNewBufSzPtr )
; real_T Assignment_1_17052017_acc_rt_TDelayInterpolate ( real_T tMinusDelay
, real_T tStart , real_T * tBuf , real_T * uBuf , int_T bufSz , int_T *
lastIdx , int_T oldestIdx , int_T newIdx , real_T initOutput , boolean_T
discrete , boolean_T minorStepAndTAtLastMajorOutput ) ; extern real_T
look1_pbinlxpw ( real_T u0 , const real_T bp0 [ ] , const real_T table [ ] ,
uint32_T prevIndex [ ] , uint32_T maxIndex ) ; extern real_T look1_binlxpw (
real_T u0 , const real_T bp0 [ ] , const real_T table [ ] , uint32_T maxIndex
) ; void Assignment_1_17052017_Saturation_Init ( SimStruct * S ,
DW_Saturation_Assignment_1_17052017_T * localDW ,
P_Saturation_Assignment_1_17052017_T * localP ) ; void
Assignment_1_17052017_Saturation_Update ( SimStruct * S , real_T rtu_Enable ,
B_Saturation_Assignment_1_17052017_T * localB ,
DW_Saturation_Assignment_1_17052017_T * localDW ) ; void
Assignment_1_17052017_Saturation ( SimStruct * S , real_T rtu_Enable , const
real_T rtu_phi [ 5 ] , B_Saturation_Assignment_1_17052017_T * localB ,
DW_Saturation_Assignment_1_17052017_T * localDW ,
P_Saturation_Assignment_1_17052017_T * localP ) ; void
Assignment_1_17052017_SaturationTID5 ( SimStruct * S ,
B_Saturation_Assignment_1_17052017_T * localB ,
P_Saturation_Assignment_1_17052017_T * localP ) ; void
Assignment_1_17052017_IfActionSubsystem2_Enable ( SimStruct * S ) ; void
Assignment_1_17052017_IfActionSubsystem2_Disable ( SimStruct * S ) ; void
Assignment_1_17052017_IfActionSubsystem2 ( SimStruct * S , real_T rtu_In1 ,
real_T * rty_Out1 ) ; void Assignment_1_17052017_NegSeqComputation_Disable (
SimStruct * S , DW_NegSeqComputation_Assignment_1_17052017_T * localDW ) ;
void Assignment_1_17052017_NegSeqComputation ( SimStruct * S , real_T
rtu_Enable , creal_T rtu_In , creal_T rtu_In_h , creal_T rtu_In_l ,
B_NegSeqComputation_Assignment_1_17052017_T * localB ,
DW_NegSeqComputation_Assignment_1_17052017_T * localDW ,
P_NegSeqComputation_Assignment_1_17052017_T * localP ) ; void
Assignment_1_17052017_ZeroSeqComputation_Disable ( SimStruct * S ,
DW_ZeroSeqComputation_Assignment_1_17052017_T * localDW ) ; void
Assignment_1_17052017_ZeroSeqComputation ( SimStruct * S , real_T rtu_Enable
, creal_T rtu_In , creal_T rtu_In_f , creal_T rtu_In_o ,
B_ZeroSeqComputation_Assignment_1_17052017_T * localB ,
DW_ZeroSeqComputation_Assignment_1_17052017_T * localDW ,
P_ZeroSeqComputation_Assignment_1_17052017_T * localP ) ;
#endif
