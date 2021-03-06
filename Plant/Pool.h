#ifndef PoolH
#define PoolH
#include <string>
#include "Biomass.h"
#include "StageBasedInterpolationFunction.h"

class ScienceAPI;
class Delta;
class Pool : public Biomass
   {
   public:
      Pool(plantInterface& plant, ScienceAPI& scienceAPI, const std::string& Name, const std::string& PartName,
           bool DoRegistrations = true);
      void Init ();

      StageBasedInterpolationFunction DigestibilityMax;
      StageBasedInterpolationFunction DigestibilityAvg;
      StageBasedInterpolationFunction DigestibilityMin;

      virtual Biomass& operator = (const Biomass& Pool2);

   protected:
      std::string PartName;
      std::string Name;
      ScienceAPI& scienceAPI;
      plantInterface& Plant;

      void DoRegistrations();
      virtual void CheckBounds();
   };

#endif
