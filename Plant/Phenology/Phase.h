#ifndef PhaseH
#define PhaseH

#include <string>

// Terminology:
// A "stage" is a point in time.
// A "phase" is the period between two stages.
// A "composite phase" is a group of one or more phases.

/////////////////////////////////////////////////////////////////////////////////////////////
// A phenological phase.
class Phase
   {
   protected:
     std::string myName;   // The name of the "stage" that the phase starts from.
     float tt,             // Thermal time spent in this phase
           target,         // Target time we want to spend here
           days;           // Number of days spent in this phase.
     float tt_after;
     float days_after;
     ScienceAPI& scienceAPI;
     plantInterface& plant;
     int DaysFromSowingToEndOfPhase;
     std::string EndStageName;
     interpolationFunction y_tt;
   public:

     Phase(ScienceAPI& api, plantInterface& p, const std::string& StartStageName, const std::string& EndStageName);
     virtual ~Phase() {};

     virtual void process() {};
     virtual void OnSow(float sowing_depth) {};

     virtual void calcPhaseDevelopment(int das,
                                       float& dlt_tt_phenol, float& phase_devel);

     void  add(float dlt_days, float dlt_tt) {days += dlt_days; tt += dlt_tt;};
     void  addToAfter(float dlt_days, float dlt_tt) {days_after += dlt_days; tt_after += dlt_tt;}
     void  setTT(float value)     {tt = value;};
     virtual float TT();
     float getTT(void) const       {return tt;};
     float getTTTarget(void) const {return target;};
     virtual void  reset(void)             {tt = days = days_after = tt_after = 0.0;}
     bool  isFirstDay(void) const  {return days <= 1.0;}
     float daysInPhase() {return days;}
     string name(void) const {return myName;};
     virtual string description(void)  {return "";};
     virtual void read();
     virtual void updateTTTargets(){};
     //virtual void onSow(protocol::SowType Sow){};
     virtual void setupTTTarget(void){};

     virtual float stress() {return 1.0;}  // no stress.
     int getDaysAfter(void) {return (int)days_after;}
     void setDaysTo(int NumDaysTo);
     std::string FullName()
        {
        string st = myName + "to" + EndStageName;
        replaceAll(st, "_", "");
        return st;
        }
   };

bool operator == (const Phase &a, const Phase &b);


#endif

