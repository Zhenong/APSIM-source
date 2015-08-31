*     ===========================================================
      subroutine Maize_Farquhar_total()
*     ===========================================================
      implicit none

*+  Sub-Program Arguments

*+  Purpose
*   Calculate daily net assimilation and dry matter production with Farquhar Model

*+  Mission Statement
*   (i)Daily assimilation is based on sunlit & shaded leaf routines
*   (ii)net assimilation minus maintenance and growth respiration is daily dm
	  
*+  Changes
*   First coded by ZN-Jin, following J.M. Chen et al(1999) Ecological Modelling

*+  Constant Values
      character  my_name*(*)           ! name of procedure
      parameter (my_name = 'Maize_Farquhar_total')
	  
*- Implementation Section ----------------------------------
	  
      call crop_radn_2leaf(
     :        g%lai
     :        g%radn
     :        g%radn_int
     :        g%dayl
     :        c%EPAR
     :        g%lai_sun
     :        g%lai_shade
     :        g%ppfd_sun_noon
     :        g%ppfd_shade_noon
     :     )	  

	! calculate sunlit leaf assimilation
      call stomatal_conductance(
     :        c%stomatal_gmax
     :        g%ppfd_sun_noon
     :        c%ppfd_coef
     :        g%maxt
     :        g%mint
     :        c%stomatal_optT
     :        c%stomatal_limT
     :        g%vpd
     :        c%vpd_close
     :        c%vpd_open
     :        g%stomatal_gn
     :        g%stomatal_gmin
     :     )

      call Farquhar_Photosynthesis(
     :        g%maxt
     :        g%mint
     :        g%nfact_photo
     :        g%ppfd_sun_noon
     :        g%stomatal_gn
     :        g%stomatal_gmin
     :        g%A_sun
     :     )	
	  
	 
	 ! calculate shaded leaf assimilation
      call stomatal_conductance(
     :        c%stomatal_gmax
     :        g%ppfd_shade_noon
     :        c%ppfd_coef
     :        g%maxt
     :        g%mint
     :        c%stomatal_optT
     :        c%stomatal_limT
     :        g%vpd
     :        c%vpd_close
     :        c%vpd_open
     :        g%stomatal_gn
     :        g%stomatal_gmin
     :     )		 

      call Farquhar_Photosynthesis(
     :        g%maxt
     :        g%mint
     :        g%nfact_photo
     :        g%ppfd_shade_noon
     :        g%stomatal_gn
     :        g%stomatal_gmin
     :        g%A_shade
     :     )	

	 
    !daily total assimilation(i.e. GPP)	 
      g%dlt_gpp = g%A_sun*g%lai_sun + g%A_shade*g%lai_shade	 
	
	! minus maintenance respiration & growth respiration
      call C4_maintRespiration(
     :        g%dm_green,
     :        g%maxt,
     :        g%mint,
     :        c%Rm_Q10,
     :        c%Rm_coef,
     :        g%dlt_Rm
	  )

      call C4_growthRespiration(
     :        g%dlt_gpp,
     :        c%Rg_coef,
     :        g%dlt_Rg

    !daily net assimilation(i.e. dry matter production)	  
      g%dlt_dm_light = g%dlt_gpp - g%dlt_Rm - g%dlt_Rg
	
      return
      end subroutine



*     ===========================================================
      subroutine Farquhar_Photosynthesis(
     :            g_maxt,
     :            g_mint,
     :            leaf_N_stress,
     :            PPFD,
     :            g_n,
     :            g_min,
     :            Anet	 
     :            )
*     ===========================================================
      implicit none

*+  Sub-Program Arguments
      real      g_maxt            !daily maximum air temperature(oC)
      real      g_mint            !daily minimum air temperature(oC)
      real      leaf_N_stress     !stress due to leaf N content
      real      PPFD              !PAR flux density, per unit projected LAI(umol photons/m2/s)
      real      g_n               !stomatal conductance at noon(umol CO2/m2/s/Pa)
      real      g_min             !stomatal conductance at minimum(umol CO2/m2/s/Pa)
      real      Anet              !(Output) daily net assimilation per leaf area(gC m-2)
	  
*+  Purpose


*+  Mission Statement
 
	  
*+  Changes
*   First coded by ZN-Jin, following J.M. Chen et al(1999) Ecological Modelling
*   Currently Jmax is assumed proportional to Vmax; in future add Temp denpendency function

*+  Constant Values
      character  my_name*(*)           ! name of procedure
      parameter (my_name = 'Farquhar_Photosynthesis')
	  
*+  Local Variables
      real    conc_O2 = 21000.0    !oxygen concentration in the atmosphere(Pa)
      real    Vmax25 = 50.0        !Vmax at 25 oC  
      real    R_gas = 8.3143       !the molar gas constant(m3 Pa mol-1 K-1)  
      real    maxN = 1.5           !maximum leaf nitrogen content(%) 	  
      real    Ca = 40.0            !Atmospheric CO2 concentration(Pa)
	  
      real    meant, stdT          !mean daily air temperature(oC)	  
      real    gamma                !CO2 compensation point without dark respiration                
      real    Kc, Ko, K
      real    f_Temp, f_Nitrogen
      real    Vmax, Jmax, J
      real    scatt = 0.15         !leaf reflectance + transmittance
      real    f_spec = 0.15        !spectral correction
      real    theta = 0.7          !curvature of leaf response of electron transport to irradiance
      real    cf_b, cf_c           !coefficients for electron transport equation
      real    I2                   !useful light absorbed by PSII(umol photons/m2/s) 
      real    Rd                   !!daytime leaf dark respiration(umol m-2 s-1)
      real    a, b, c, d, e        !intermediate variable when calculate assimilation
      real    Ac, Aj
	  
*- Implementation Section ----------------------------------
	  
	  !Standarized air temperature (relative to 25 oC)
      meant = (maxt + mint)/2.0
      stdT = (meant-25.0)/10.0
	  
      !calc gamma, CO2 compensation point without dark respiration
      gamma = 1.92 * 10**(-4) * conc_O2 * 1.75**stdT      
	  
      !calc K, enzyme kinetics
	  Kc = 30.0 * 2.1**stdT    !Michaelis-Menten constant for CO2(Pa)
	  Ko = 30000 * 1.2**stdT   !Michaelis-Menten constant for O2(Pa)
      K = Kc * (1 + conc_O2/Ko)

      !calc Vmax, maximum carboxylation rate(mu mol m-2 s-1)    
      f_Temp = (1.0 + exp((-220000.0+710.0*(meanT+273.0))/
     :      (R_gas*(meanT+273.0)))) ** (-1)
	 
      f_Nitrogen = bound(leaf_N_stress, 0.8, 1.0)  
      
      Vmax = Vmax25 * 2.4**stdT * f_Temp * f_Nitrogen

      !calc J, electron transport rate(mu mol m-2 s-1)
      !essentially de Pury & Farquhar 1997, also see MAIZSIM CGag_exchange::Photosynthesis 
      Jmax = 2.1 * Vmax                        !light-saturated rate of J
      I2 = PPFD * (1 - scatt) * (1 - f)/2.0    !useful light absorbed by PSII
      cf_b = -(I2 + Jmax)
      cf_c = I2 * Jmax
      if((cf_b**2 - 4.0*theta*cf_c).ge.0) then
          J = (-cf_b -sqrt(cf_b**2 - 4.0*theta*cf_c))/(2.0*theta)
      else                                     !in case no real solution to equation
          J = Jmax * PPFD / (PPFD + 2.1*Jmax)
      endif  		 
	  
      !Instantaneous leaf-level assimilation	  
      !Ac = Vmax * (Ci - gamma)/(Ci + K)	  
	  !Aj = J * (Ci - gamma)/(4.5*Ci + 10.5*gamma)
      Rd = 0.015 * Vmax                        !daytime leaf dark respiration	  
      !Anet = min(Ac, Aj) - Rd  

	  
!Leaf level daily intergration
      !For Ac, Eq-16 in Chen et al. (1999)
      a = (K + Ca)**2.0
      b = 2.0*(2.0*gamma + K - Ca)*Vmax + 2.0*(Ca + K)*Rd
      c = (Vmax - Rd)**2.0
      d = (a*g_n**2.0 + b*g_n + c)**0.5
      e = (a*g_min**2.0 + b*g_min + c)**0.5
	  
      Ac = 1.27/(2.0*(g_n - g_min))*(0.5*sqrt(a)*(g_n**2 - g_min**2.0) +
     :        sqrt(c)*(g_n - g_min**2.0) - d*(2.0*a*g_n + b)/(4.0*a) +
     :        sqrt(e)*(2.0*a*g_min + b)/(4.0*a) +
     :        (b**2.0 - 4.0*a*c)/(8.0*a**1.5) *
     :        log((2*a*g_n+b+2*sqrt(a)*d)/(2*a*g_min+b+2*sqrt(a)*d))) 

      !For Aj
      a = (2.3*gamma + Ca)**2.0
      b = 0.4*(4.3*gamma - Ca)*J + 2*(Ca + 2.3*gamma)*Rd
      c = (0.2*J - Rd)**2
      d = (a*g_n**2 + b*g_n + c)**0.5
      e = (a*g_min**2 + b*g_min + c)**0.5
	  
      Aj = 1.27/(2.0*(g_n - g_min))*(0.5*sqrt(a)*(g_n**2 - g_min**2.0) +
     :        sqrt(c)*(g_n - g_min**2.0) - d*(2.0*a*g_n + b)/(4.0*a) +
     :        sqrt(e)*(2.0*a*g_min + b)/(4.0*a) +
     :        (b**2.0 - 4.0*a*c)/(8.0*a**1.5) *
     :        log((2*a*g_n+b+2*sqrt(a)*d)/(2*a*g_min+b+2*sqrt(a)*d))) 	 

      !net daily assimilation, Rd was already deducted 
      Anet = min(Ac, Aj)
	  
      return
      end subroutine
	  
	  

*     ===========================================================
      subroutine crop_radn_2leaf(
     :            g_lai,
     :            radn,
     :            radn_int,
     :            g_dayl,
     :            c_EPAR,
     :            g_lai_sun,
     :            g_lai_shade,
     :            g_ppfd_sun_noon,
     :            g_ppfd_shade_noon 
     :            )
*     ===========================================================
      implicit none

*+  Sub-Program Arguments
      real      g_lai             !projected total LAI
      real      radn              !downward shortwave irradiance(MJ/m2/day)
      real      radn_int          !leaf intercepted radiation based on light extinction(MJ/m2/day)
      real      g_dayl            !day length(hr)
      real      c_EPAR            !PAR photon energy ratio = 4.55(umol/J)
      real      g_lai_sun         !(Output) projected sunlit leaf area index(m2/m2)
      real      g_lai_shade       !(Output) projected shaded leaf area index(m2/m2) 
      real      g_ppfd_sun_noon   !(Output) ppfd for sunlit leaf at noon(umol/m2/s)
      real      g_ppfd_shade_noon !(Output) ppfd for shaded leaf at noon(umol/m2/s)
	  
*+  Purpose
*   Seek radiance absorbed by sunlit & shaded leaf

*+  Mission Statement
*   (i) simulate sunlit & shaded LAI
*   (ii) distribute downward shortwave radiance flux(MJ/m2/day) to sunlit/shaded leaf component  
	  
*+  Changes
*   First coded by ZN-Jin, following Biome-BGC 4.2 ("radtrans.c")

*+  Constant Values
      character  my_name*(*)           ! name of procedure
      parameter (my_name = 'crop_radn_2leaf')
	  
*+  Local Variables
      real      ext_coef      !canopy light extinction coefficient, possibly already given in APSIM
      real      PARabs        !absorbed PAR by whole canopy(W/m2)
      real      PARabs_sun    !absorbed PAR by sunlit leaf(W/m2)
      real      PARabs_shade  !absorbed PAR by shaded leaf(W/m2)
      real      g_ppfd_per_sun    !daily average ppfd for sunlit leaf(umol/m2/s)
      real      g_ppfd_per_shade  !daily average ppfd for shaded leaf(umol/m2/s)
	  
*- Implementation Section ----------------------------------

      g_lai_sun = 1 - exp(- g_lai)    !sunlit LAI will never exceed 1.0
      g_lai_shade = g_lai - g_lai_sun	  

	!absorbed total solar irradiance
      ext_coef = log(1 - radn_int/radn) / (-g_lai)
      PARabs = radn_int / 0.0864    !convert MJ/m2/day to W/m2
	
    !decompose to sunlit & shaded components	
      PARabs_sun = ext_coef * (radn/0.0864) * g_lai_sun    !0.0864 for unit conversion
      PARabs_shade = l_bound(PARabs - PARabs_sun, 0.0)

    !convert to PPFD: assumes an average energy for PAR photon (EPAR, umol/J)
    !first calculate per unit area based radiation
    !unit conversion: W/m2 --> umol/m2/s	
      g_ppfd_per_sun = (PARabs_sun / g_lai_sun) * c_EPAR
      g_ppfd_per_shade = (PARabs_shade / g_lai_shade) * c_EPAR
	  
    !PPFD at noon  
      g_ppfd_sun_noon = g_ppfd_per_sun * (3.1415/2.0) / (g_dayl/24.0)
      g_ppfd_shade_noon = g_ppfd_per_shade * (3.1415/2.0) / (g_dayl/24.0)
	  
      return
      end subroutine

	  
	  
*     ===========================================================
      subroutine stomatal_conductance(
     :            c_gmax,
     :            ppfd,
     :            ppfd_coef,
     :            g_maxt,
     :            g_mint,
     :            c_Topt,
     :            c_Tlim,
     :            g_vpd,
     :            c_vpd_close,
     :            c_vpd_open,
     :            g_n,
     :            g_min
     :            )
*     ===========================================================
      implicit none

*+  Sub-Program Arguments
      real      c_gmax       !maximum stomatal conductance(m s-1)?
      real      ppfd         !photosynthetic photo flux density(umol/m2/s)
      real      ppfd_coef    !coef in a relationship between gs and ppfd(umol/m2/s)
      real      g_maxt       !daily maximum temperature(oC)
      real      g_mint       !daily minimum temperature(oC)
      real      c_Topt       !optimum temperature(oC)
      real      c_Tlim       !limiting temperature beyond which no conductance(oC)
      real      g_vpd        !daily average VPD(kPa)
      real      c_vpd_close  !VPD at stomatal closure(kPa)
      real      c_vpd_open   !VPD at stomatal full opening(kPa)
      real      g_n          !(Output) stomatal conductance at noon(umol m-2 s-1 Pa-1)
      real      g_min        !(Output) minimum stomatal conductance(umol m-2 s-1 Pa-1)
	  
*+  Purpose
*   calculate stomatal conductance(g) at noon(g_n) & g_min

*+  Mission Statement
 
	  
*+  Changes
*   First coded by ZN-Jin, following J.M. Chen et al(1999) Ecological Modelling

*+  Constant Values
      character  my_name*(*)           ! name of procedure
      parameter (my_name = 'stomatal_conductance')
	  
*+  Local Variables
      real      meant     !daily mean temperature(oC) 
      real      f_ppfd    !multiplier of radiation stress
      real      f_temp    !multiplier of temperature stress
      real      f_vpd     !multiplier of VPD stress
      real      R_gas = 8.3143       !the molar gas constant(m3 Pa mol-1 K-1) 
  	  
*- Implementation Section ----------------------------------
      
	  !radiation stress
      f_ppfd = ppfd * ppfd_coef / (1 + ppfd*ppfd_coef)

      !temperature stress	  
      meant = (g_maxt + g_mint)/2.0
      if(meant.le.c_Topt .and. meant.lt.1.0) then
          f_Temp = log(meant)/log(c_Topt)
      elseif(meant.gt.c_Topt) then
          f_Temp = cos(3.1415/2.0*(meant - c_Topt)/(c_Tlim - c_Topt))
      else
          f_Temp = 0.0  
      endif

      !VPD stress
      if(g_vpd.le.c_vpd_open) then
          f_vpd = 1.0
      elseif(g_vpd.gt.c_vpd_open .and. g_vpd.lt.c_vpd_close) then
          f_vpd = (c_vpd_close - g_vpd)/(c_vpd_close - c_vpd_open)
      else
          f_vpd = 0.0  
      endif	  

      !stomatal conductance at noon and minimum	  
      g_n = c_gmax * f_ppfd * f_Temp * f_vpd
      g_min = 0.0        !prescribed
	  
	  !Be careful, now conductance unit in (m s-1)
	  !Need convert to umol m-2 s-1
      g_n = g_n * 1e6 / (R_gas*(meant + 273.0)) 
      g_min = g_min * 1e6 / (R_gas*(meant + 273.0))
	  
      return
      end subroutine

	  
	  
*     ===========================================================
      subroutine C4_maintRespiration(
     :            Biomass,
     :            g_maxt,
     :            g_mint,
     :            c_Rm_Q10,
     :            c_Rm_coef,
     :            Rm	 
     :            )
*     ===========================================================
      implicit none

*+  Sub-Program Arguments
      real    Biomass        !current day total biomass(i.e. dry matter)(gC m-2)
      real    g_maxt         !daily maximum air temperature(oC)
      real    g_mint         !daily minimum air temperature(oC)
      real    c_Rm_Q10       !Q10 value for maintenance respiration function
      real    c_Rm_coef      !coefficient for maintenance respiration(0.002 g g-1 day-1 at 20oC)
      real    Rm             !daily maintenance respiration(gC m-2 day-1)
	  
*+  Purpose
*   Calculate daily maintenance respiration

*+  Mission Statement
*   Maintenance respiration is a function of current biomass and temperature
	  
*+  Changes
*   First coded by ZN-Jin, calculation rounte from MAIZSIM, parameters from JM Chen(1999) 

*+  Constant Values
      character  my_name*(*)           ! name of procedure
      parameter (my_name = 'C4_maintRespiration')
	  
*+  Local Variables
      real    baseT = 20.0        !base temperature for Q10 function(oC)
      real    meant               !daily mean air temperature(oC)
	  
*- Implementation Section ----------------------------------
    
      meant = (maxt + mint)/2.0
	  
      Rm = Biomass * c_Rm_coef * c_Rm_Q10**((meant - baseT)/10.0)
	  
      return
      end subroutine


	  
*     ===========================================================
      subroutine C4_growthRespiration(
     :            GPP,
     :            c_Rg_coef,
     :            Rg	 
     :            )
*     ===========================================================
      implicit none

*+  Sub-Program Arguments
      real      GPP         !Gross Primary Production(gC m-2 day-1)
      real      c_Rg_coef   !overall growth respiration coefficient(0-1)	  
      real      Rg          !Growth respiration(gC m-2 day-1)
	 
*+  Purpose
*   Calculate growth respiration

*+  Mission Statement
*   Growth respiration is proportional to daily GPP
*   We assume single coefficient for overall Rg, but can be potentially
*   further divided into different plant parts
	  
*+  Changes
*   First coded by ZN-Jin, following MAIZSIM 

*+  Constant Values
      character  my_name*(*)           ! name of procedure
      parameter (my_name = 'C4_growthRespiration')
	  
*+  Local Variables
	  
*- Implementation Section ----------------------------------

      Rg = GPP * c_Rg_coef
	  
      return
      end subroutine
