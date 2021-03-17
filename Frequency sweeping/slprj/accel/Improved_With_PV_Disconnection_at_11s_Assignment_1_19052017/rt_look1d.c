

    /*
  * rt_look1d.c
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



    

    #include "rt_look1d.h"



  

  

  

  

  

  

  

  

  

  

  

  

  

  

  

    
  /* 1D lookup routine for data type of real_T. */
     real_T rt_Lookup(const real_T *x, int_T xlen, real_T u, const real_T *y)
   {
     int_T idx = rt_GetLookupIndex(x, xlen, u);
     real_T num = y[idx+1] - y[idx];
     real_T den = x[idx+1] - x[idx];
     
     /* Due to the way the binary search is implemented
     in rt_look.c (rt_GetLookupIndex), den cannot be
     0.  Equivalently, m cannot be inf or nan. */
     
     real_T m = num/den;
     
     return (y[idx] + (m * (u - x[idx])));
   }
     
  


  

  

  

  
