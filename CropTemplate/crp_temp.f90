      Module crp_tempModule

      contains

!     ===========================================================
      subroutine crop_temperature_stress_photo (num_ave_temp          &
                ,x_ave_temp,y_stress_photo, maxt, mint,          &
                temp_stress_photo)
!     ===========================================================

!      dll_export crop_temperature_stress_photo
      use scienceModule
      use dataModule
      use errorModule
      implicit none

!+  Sub-Program Arguments
      INTEGER num_ave_temp        ! (INPUT)  size_of critical temperature table
      REAL    x_ave_temp(*)       ! (INPUT)  critical temperatures for photosynthesis (oC)
      REAL    y_stress_photo(*)   ! (INPUT)  Factors for critical temperatures (0-1)
      REAL    maxt                ! (INPUT)  maximum air temperature (oC)
      REAL    mint                ! (INPUT)  minimum air temperature (oC)
      REAL    temp_stress_photo   ! (OUTPUT) photosynthetic reduction factor for
                                  ! temperature stress (0-1)

!+  Purpose
!       photosynthetic reduction factor for temperature stress (0-1)

!+  Mission Statement
!   Calculate the temperature factor for photosynthesis

!+  Changes
!       090994 jngh specified and programmed
!       970216 slw generalised

!+  Constant Values
      character  my_name*(*)      ! name of procedure
      parameter (my_name = 'crop_temperature_stress_photo')

!+  Local Variables
      real       ave_temp         ! mean temperature for the day (oC)

!- Implementation Section ----------------------------------

      call push_routine (my_name)

         ! now get the temperature stress factor that reduces
         ! photosynthesis (0-1)
         ave_temp = (maxt + mint) /2.0
         temp_stress_photo = linear_interp_real (ave_temp          &
                          , x_ave_temp, y_stress_photo          &
                          , num_ave_temp)
         temp_stress_photo = bound (temp_stress_photo, 0.0, 1.0)


      call pop_routine (my_name)
      return
      end subroutine

	  
!     ===========================================================
      subroutine crop_canopy_temperature (maxt, mint, rnet,         &
                lai, es, transpiration, height,                       & 
                canopy_temp)
!     ===========================================================

!      dll_export crop_temperature_stress_photo
      use scienceModule
      use dataModule
      use errorModule
      implicit none

!+  Sub-Program Arguments
      REAL    maxt                ! (INPUT)  maximum air temperature (oC)
      REAL    mint                ! (INPUT)  minimum air temperature (oC)
      REAL    rnet                ! (INPUT)  net daily radiation (MJ/m2)
      REAL    lai                 ! (INPUT)  leaf area index of previous day (m2/m2)
      REAL    es                  ! (INPUT)  daily soil evaporation (mm)
      REAL    transpiration       ! (INPUT)  daily plant transpiration (mm)
      REAL    height              ! (INPUT)  plant height (mm), need convert to (m)
      REAL    canopy_temp         ! (OUTPUT) canopy temperature based on STICS empirical method

!+  Purpose
!       photosynthetic reduction factor for temperature stress (0-1)
!       canopy temperature according to STICS empirical approach instead of air temperature is used

!+  Mission Statement
!   Calculate the temperature factor for photosynthesis

!+  Changes
!   created by ZN-J

!+  Constant Values
      character  my_name*(*)      ! name of procedure
      parameter (my_name = 'crop_temperature_stress_photo')

!+  Local Variables
      real       ave_temp         ! mean temperature for the day (oC)
      real       tcult_max		  ! daily maximum canopy temperature (oC)
      real       tcult_min		  ! daily maximum canopy temperature (oC)
      real       Z0		          ! empirical fraction of canopy height (m)
	  
!- Implementation Section ----------------------------------

      call push_routine (my_name)
         ! calculate maximum canopy temperature
         if (lai .ge. 1.5) then                   ! when plant is alive and canopy is dense
             Z0 = max(0.13*height/1000, 0.001)	 
             tcult_max = maxt + (rnet/2.46 - es - transpiration -1.27) &
                                /(1.68/log(1/Z0))
             tcult_min = mint
             ave_temp = (tcult_max+ tcult_min) /2.0
         else                                     ! when plant is out or canopy is not dense
		     ave_temp = (maxt+ mint) /2.0
         endif
		 
         canopy_temp = ave_temp
		 
      call pop_routine (my_name)
      return
      end subroutine

	  
	  

!     ===========================================================
      subroutine crop_temperature_stress_swat (baseT          &
                ,optT, maxt, mint, temp_stress_photo)
!     ===========================================================

!      dll_export crop_temperature_stress_photo
      use scienceModule
      use dataModule
      use errorModule
      implicit none

!+  Sub-Program Arguments
      REAL    baseT       		  ! (INPUT)  base temperatures for photosynthesis (oC)
      REAL    optT   			  ! (INPUT)  optimum temperatures for photosynthesis (oC)
      REAL    maxt                ! (INPUT)  maximum air temperature (oC)
      REAL    mint                ! (INPUT)  minimum air temperature (oC)
      REAL    temp_stress_photo   ! (OUTPUT) photosynthetic reduction factor for
                                  ! temperature stress (0-1)

!+  Purpose
!        photosynthetic reduction factor for temperature stress (0-1)

!+  Mission Statement
!   Calculate the temperature factor for photosynthesis with SWAT method

!+  Changes
!	First added by ZN-J

!+  Constant Values
      character  my_name*(*)      ! name of procedure
      parameter (my_name = 'crop_temperature_stress_photo')

!+  Local Variables
      real       ave_temp         ! mean temperature for the day (oC)
      real       tstrs            ! temporary variable for temperature stress (0-1)

!- Implementation Section ----------------------------------

      call push_routine (my_name)

         ! now get the temperature stress factor that reduces photosynthesis (0-1)
         ave_temp = (maxt + mint) /2.0
 
         if (ave_temp.le.baseT) then
                 tstrs = 1
 
         elseif (ave_temp.gt.baseT .and. ave_temp.le.optT) then
             tstrs = 1-exp((-0.1054*(optT - ave_temp)**2.0)/   &
                          ((ave_temp - baseT)*2.0)) 
         
         elseif (ave_temp.gt.optT .and. ave_temp.le.(2*optT-baseT))then
             tstrs = 1-exp((-0.1054*(optT - ave_temp)**2.0)/   &
                          ((2.0*optT - ave_temp - baseT)*2.0)) 
         else
             tstrs = 1
         endif
  
         temp_stress_photo = bound (tstrs, 0.0, 1.0)


      call pop_routine (my_name)
      return
      end subroutine

  
!     ===========================================================
      real function crop_htstress_HSA(g_maxt, g_mint,            &
                                   c_HSA_crT, c_HSA_limT)
!     ===========================================================	  
      use dataModule                 ! for bound check function
      implicit none

!+ Sub-Program Arguments
      real       g_maxt              ! (INPUT) maximum daily air temperature
      real       g_mint              ! (INPUT) minimum daily air temperature
      real       c_HSA_crT           ! (INPUT) critical temp above which grain-set decline
      real       c_HSA_limT          ! (INPUT) limiting temp above which grain-set stop

!+ Purpose
!  calculate heat stress on HI with PEGASUS method

!+ Mission Statement
!  calculate a linearly decease scalar to HI between c_HSA_crT and c_HSA_limT
!  algorithm from PEGASUS model, but first introduced by Challinor (2004) in GLAM model

!+ Changes
!  First added by ZN-J

!+ Constant Values
      character  myname*(*)            ! procedure name
      parameter (myname = 'crop_htstress_HSA')

!+ Local Variables
      real      am_temp                ! daily effective temperature (8 am to 2 pm)
									   ! am_temp = (Tmax + Tmean)/2
      real      ht_stress              ! daily heat stress
	  
!- Implementation Section ----------------------------------
      am_temp = 0.75*g_maxt + 0.25*g_mint
	  
	  ht_stress = 1-(am_temp - c_HSA_crT )/(c_HSA_limT - c_HSA_crT)

      crop_htstress_HSA = bound(ht_stress, 0.0, 1.0)      

      return
      end function


      end module crp_tempModule
