

    /*
  * rt_look.h
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


  #ifndef RTW_HEADER_rt_look_h_
  #define RTW_HEADER_rt_look_h_
  

    

    #include "rtwtypes.h"
  
  #ifdef DOINTERPSEARCH
  #include <float.h>
  #endif


  

  

  

  

  

  

  

  

  

  

  

  

  

  

    
    
  #ifndef INTERP
  # define INTERP(x,x1,x2,y1,y2) ( (y1)+(((y2) - (y1))/((x2) - (x1)))*((x)-(x1)) )
  #endif

  #ifndef ZEROTECHNIQUE
  #define ZEROTECHNIQUE
  typedef enum {
    NORMAL_INTERP,
    AVERAGE_VALUE,
    MIDDLE_VALUE
  } ZeroTechnique;
  #endif



    extern int_T rt_GetLookupIndex(const real_T *x, int_T xlen, real_T u) ;
  
  


  

  

  

  

  

    #endif /* RTW_HEADER_rt_look_h_ */
