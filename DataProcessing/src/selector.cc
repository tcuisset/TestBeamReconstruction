#include "TestBeamReconstruction/DataProcessing/interface/selector.h"
#include "TCanvas.h"

Selector::Selector(const std::string& in_file_path, const std::string& out_file_path,
		   const std::string& datatype_, const std::string& showertype_, const int& beam_energy_,
		   std::optional<std::string> in_tree_name, std::optional<std::string> in_tree_name_friend, std::optional<std::string> out_tree_name)
{
  sanity_checks(in_file_path);
  this->indata_.file_path = in_file_path;
  this->outdata_.file_path = out_file_path;
  this->beam_energy = beam_energy_;

  if(datatype_.find("data") != std::string::npos)
    this->datatype = DATATYPE::DATA;
  else if(datatype_.find("sim") != std::string::npos)
    this->datatype = DATATYPE::MC;
  else
    throw std::invalid_argument("Wrong data type: " + datatype_);
    
  if(showertype_ == "em")
    this->showertype = SHOWERTYPE::EM;
  else if(showertype_ == "had")
    this->showertype = SHOWERTYPE::HAD;
  
  if(in_tree_name != std::nullopt)
    this->indata_.tree_name = in_tree_name.value();
  if(in_tree_name_friend != std::nullopt)
    this->indata_.tree_name_friend = in_tree_name_friend.value();
  if(out_tree_name != std::nullopt)
    this->indata_.tree_name = out_tree_name.value();

  //establish which data columns will be saved
  this->savedcols_ = {"event", "run", "NRechits", new_detid_, "ce_clean_x_unshifted", "ce_clean_y_unshifted", new_z_, new_layer_, 
    new_en_, new_en_MeV_, new_ahc_en_MeV_, "beamEnergy",
    "impactX_unshifted", "impactY_unshifted", "DWC_b_x", "DWC_b_y", "DWC_trackChi2_X", "DWC_trackChi2_Y"};
  for(unsigned i=1; i<=detectorConstants::totalnlayers; ++i) {
    impactXcols_.push_back("myFriend.impactX_HGCal_layer_" + std::to_string(i));
    impactcols_.push_back("myFriend.impactX_HGCal_layer_" + std::to_string(i));
    impactYcols_.push_back("myFriend.impactY_HGCal_layer_" + std::to_string(i));
    impactcols_.push_back("myFriend.impactY_HGCal_layer_" + std::to_string(i));
  }
  

  //load external data
  this->load_noise_values();
  this->load_shift_values();
}

Selector::~Selector()
{
}

void Selector::load_noise_values()
{
  std::string rel( getenv("CMSSW_BASE") );
  std::string filename = rel + "/src/TestBeamReconstruction/DataProcessing/Noise_Map.txt";
  std::ifstream in(filename);
  if(!in.is_open())
    throw std::invalid_argument("File " + filename + " could not be opened.");

  std::pair<int,int> mod_chip;
  std::pair<std::pair<int,int>, float> temp;
  unsigned layer, mod_id, mod_pos, chip;
  float noise;

  for(unsigned i=0; i<2; ++i) //lines to skip
    in.ignore(std::numeric_limits<unsigned>::max(), '\n');
  while(in >> layer >> mod_id >> mod_pos >> chip >> noise)
    {
      mod_chip = std::make_pair(mod_id, chip);
      temp = std::make_pair(mod_chip, noise);
      this->noise_map_.insert(temp);
    }
}

void Selector::load_shift_values()
{
  std::string rel( getenv("CMSSW_BASE") );
  std::string filename = rel + "/src/TestBeamReconstruction/DataProcessing/Impact_Shifts.txt";
  std::ifstream in(filename);
  if(!in.is_open())
    throw std::invalid_argument("File " + filename + " could not be opened.");

  unsigned layer;
  float shiftx, shifty;
  for(unsigned i=0; i<2; ++i) //lines to skip
    in.ignore(std::numeric_limits<unsigned>::max(), '\n');
  for(unsigned i=0; i<detectorConstants::totalnlayers; ++i) {
    in >> layer >> shiftx >> shifty;
    assert(layer == i+1);
    this->shifts_map_.emplace_back( shiftx, shifty );
  }
}

/**
 * Lookup noise map and find whether we are in nsigma of noise of the chip (3 sigma for em, 4 for hadronic shower)
*/
bool Selector::reject_noise(const mapT& map, const unsigned& mod, const unsigned& chip, const unsigned& l, const float& amp, const bool& st)
{
  float temp_noise = -1.0;
  std::pair<unsigned, unsigned> mod_chip = std::make_pair(mod, chip);
  mapT::const_iterator it = map.find(mod_chip);
  if(it != map.end())
    temp_noise = it->second;
  else
    throw std::out_of_range("Value NOT found for Module = " + std::to_string(mod_chip.first) + " and chip = " + std::to_string(mod_chip.second));

  float sigma = 3.;
  if(st)
    sigma = 4.;
  
  return amp < sigma * temp_noise;
}

/**
 * For a given rechits, whose info is given in parameters, determine whether it passes selections
 * \param map the noiseMap (same for all hits, see \a Selector::noise_map_
*/
bool Selector::common_selection(const unsigned& layer, const float& energy, const unsigned& chip, const unsigned& channel, const unsigned& module, const float& amplitude, const bool& noise_flag, const mapT& map, const bool& showertype)
{
  if(layer<1 || layer>detectorConstants::totalnlayers)
    throw std::out_of_range("Unphysical layer number: " + std::to_string(layer));
  else if(!showertype and layer>detectorConstants::nlayers_emshowers) 
    return false; //If we want em showers (showertype=false) then only consider first 28 layers 

  if(noise_flag==1)
    return false;

  if(energy<=0.5) //MIP noise cut
    return false;
  if( (chip==0 and layer==1) or (chip==0 and channel==44) or (chip==3 and channel==22) or (chip==3 and channel==28) ) //mask noisy chips
    return false;
  if( chip==1 and layer==36 and channel==16 and module==39 )
    return false;

  if( reject_noise(map, module, chip, layer, amplitude, showertype) )
    return false;
    
  return true;
}

/**
 * cleans hits in CE-E and CE-H
 * The params are : first var is variable of interest
 * Then all columns of clean_cols ("rechit_energy", "rechit_layer", "rechit_chip", "rechit_channel", "rechit_module", "rechit_amplitudeHigh", "rechit_noise_flag", "st")
 *  with the noise map added before "st"
 * This will select from param var only the entries which passes selections, and return the subset of values
 */
template<typename T>
std::vector<T> Selector::clean_ce(const std::vector<T>& var, const std::vector<float>& en, const std::vector<unsigned>& l, const std::vector<unsigned>& chip, const std::vector<unsigned>& channel, const std::vector<unsigned>& module, const std::vector<float>& amplitude, const std::vector<bool>& noise_flag, const mapT& map, const bool& st)
{
  size_t nhits = en.size();
  std::vector<T> var_clean;
  var_clean.reserve(nhits);

  for(unsigned i=0; i<nhits; ++i) {
    if( common_selection(l[i], en[i], chip[i], channel[i], module[i], amplitude[i], noise_flag[i], map, st) )
      var_clean.push_back( var[i] );
  }
  return var_clean;  
}

//cleans hits in AHCAL
template<typename T>
std::vector<T> Selector::clean_ahc(const std::vector<T>& var, const std::vector<float>& en, const std::vector<int>& l, const bool& st)
{
  size_t nhits = var.size();
  std::vector<T> var_clean;
  var_clean.reserve(nhits);

  if(st)
    {
      for(unsigned i=0; i<nhits; ++i) {
	if(l[i] != 38) //mask AHCAL layer #38
	  var_clean.push_back(var[i]);
      }
    }
  else //no change for em showers
    return var;
  return var_clean;
}

float Selector::ahc_energy_sum(const std::vector<float>& ahc_en)
{
  return std::accumulate(ahc_en.begin(), ahc_en.end(), 0.f);
}

/** 
 * weights the energy of the CE-E and CE-H hits for a single event
 * \param en vector of ce_clean_energy for all rechits that pass selections (in MIP)
 * \param l same but layer
 * \param st showertype (false= em shower, true = hadronic shower)
 */
std::vector<float> Selector::weight_energy_ce(const std::vector<float>& en, const std::vector<unsigned>& l, const bool& st) 
{
  unsigned layer = 0;
  float energy = 0.f;
  float weight = 1.f;
  size_t nhits = en.size();
  std::vector<float> en_weighted;
  en_weighted.reserve(nhits); //maximum possible size

  for(unsigned i=0; i<nhits; ++i)
    {
      energy = en[i];
      layer = l[i];

      if( (!st and layer>0 and layer<=detectorConstants::nlayers_emshowers) //em shower case
        or (st and layer>0 and layer<=detectorConstants::totalnlayers) ) //had case
      {
        if(st and layer>detectorConstants::nlayers_emshowers) //had
          weight = detectorConstants::globalWeightCEH;
        else //em shower or had shower in CE-E
          weight = detectorConstants::dEdX.at(layer-1);
      }
      else
        continue;

      en_weighted.push_back( weight * energy ); 
    }
  return en_weighted;  
}

//weights the energy of the AHCAL hits
std::vector<float> Selector::weight_energy_ahc(const std::vector<float>& en, const bool& st) 
{
  float weight = detectorConstants::globalWeightCEH * detectorConstants::globalWeightRelative;
  size_t nhits = en.size();
  std::vector<float> en_weighted;
  en_weighted.reserve(nhits); //maximum possible size

  if(!st) // em shower
    weight = 0.f;

  for(unsigned i=0; i<nhits; ++i)
    en_weighted.push_back( weight * en[i] );
      
  return en_weighted;  
}

/**
 * \param v vector of values of columns :impactX_HGCal_layer_i, impactY_HGCal_layer_i for i in all layers
*/
bool Selector::remove_missing_dwc(const std::vector<float>& v)
{
  for(auto& x : v) {
    if(x == -999)
      return false;
  }
  return true;
}
    
void Selector::select_relevant_branches()
{
  ROOT::RDataFrame *d;
  TFile *f_had = nullptr;
  TTree *t_had1 = nullptr;
  TTree *t_had2 = nullptr;

  f_had = TFile::Open(this->indata_.file_path.c_str());
  t_had1 = static_cast<TTree*>( f_had->Get(this->indata_.tree_name.c_str()) ); ///< rechitntupler/hits TTree
  t_had2 = static_cast<TTree*>( f_had->Get(this->indata_.tree_name_friend.c_str()) ); ///< trackimpactntupler/impactPoints TTree
  t_had1->AddFriend(t_had2, "myFriend");
  d = new ROOT::RDataFrame(*t_had1);

  ROOT::Detail::RDF::ColumnNames_t clean_cols = {"rechit_energy", "rechit_layer", "rechit_chip", "rechit_channel", "rechit_module", "rechit_amplitudeHigh", "rechit_noise_flag", "st"};
  ROOT::Detail::RDF::ColumnNames_t clean_cols_detid = clean_cols; clean_cols_detid.insert( clean_cols_detid.begin(), "rechit_detid");
  ROOT::Detail::RDF::ColumnNames_t clean_cols_x = clean_cols; clean_cols_x.insert( clean_cols_x.begin(), "rechit_x");
  ROOT::Detail::RDF::ColumnNames_t clean_cols_y = clean_cols; clean_cols_y.insert( clean_cols_y.begin(), "rechit_y");
  ROOT::Detail::RDF::ColumnNames_t clean_cols_z = clean_cols; clean_cols_z.insert( clean_cols_z.begin(), "rechit_z");
  ROOT::Detail::RDF::ColumnNames_t clean_cols_layer = clean_cols; clean_cols_layer.insert( clean_cols_layer.begin(), "rechit_layer");
  ROOT::Detail::RDF::ColumnNames_t clean_cols_en = clean_cols; clean_cols_en.insert( clean_cols_en.begin(), "rechit_energy");
  
  /** wrappers required to circumvent the static-ness of the wrapped methods (required by RDataFrame) while passing a class variable (this->noise_map_)
   * \param var the variable
   * \param allOthers columns from clean_cols
  */
  auto wrapper_float = [this](const std::vector<float>& var, const std::vector<float>& en, const std::vector<unsigned>& l, const std::vector<unsigned>& chip, const std::vector<unsigned>& channel, const std::vector<unsigned>& module, const std::vector<float>& amplitude, const std::vector<bool>& noise_flag, const bool& st)
			{
			  std::vector<float> v = clean_ce<float>(var, en, l, chip, channel, module, amplitude, noise_flag, this->noise_map_, st);
			  return v;
			};
  auto wrapper_uint = [this](const std::vector<unsigned>& var, const std::vector<float>& en, const std::vector<unsigned>& l, const std::vector<unsigned>& chip, const std::vector<unsigned>& channel, const std::vector<unsigned>& module, const std::vector<float>& amplitude, const std::vector<bool>& noise_flag, const bool& st)
			{
			  std::vector<unsigned> v = clean_ce<unsigned>(var, en, l, chip, channel, module, amplitude, noise_flag, this->noise_map_, st);
			  return v;
			};


  std::string filter_str = "true", define_str = "return false;";
  if(this->showertype == SHOWERTYPE::HAD) {
    filter_str = "myFriend.dwcReferenceType>=13 && myFriend.trackChi2_X<10 && myFriend.trackChi2_Y<10";
    define_str = "return true;";
  }

  /* From ROOT documentation :
  PassAsVec is a callable generator that allows passing N variables of type T to a function as a single collection.
  PassAsVec<N, T>(func) returns a callable that takes N arguments of type T, passes them down to function func as an initializer list {t1, t2, t3,..., tN} and returns whatever f({t1, t2, t3, ..., tN}) returns.
  Note that for this to work with RDataFrame the type of all columns that the callable is applied to must be exactly T. Example usage together with RDataFrame ("varX" columns must all be float variables):
  bool myVecFunc(std::vector<float> args);
  df.Filter(PassAsVec<3, float>(myVecFunc), {"var1", "var2", "var3"});
  */
  auto partial_process = d->Filter(filter_str.c_str())
    // Remove events with missing DWC info (if any column impactX_HGCal_layer_i or  impactY_HGCal_layer_i is unset)
    .Filter(ROOT::RDF::PassAsVec<static_cast<unsigned>(2*detectorConstants::totalnlayers), float>(remove_missing_dwc),
	    impactcols_)

    .Alias("DWC_b_x", "myFriend.b_x") //Track offsets from Delay Wire Chambers : impact on EE
    .Alias("DWC_b_y", "myFriend.b_y")
    .Alias("DWC_trackChi2_X", "myFriend.trackChi2_X") //Track chisquare from DWC
    .Alias("DWC_trackChi2_Y", "myFriend.trackChi2_Y")

    // clean_cols.back() is "st" ie showertype. false -> em shower, true -> hadronic shower (initialized depending on command line argument)
    .Define(clean_cols.back(), define_str)

    //Define new columns as vector for each event, selecting only in each event rechits that pass selections
    .Define(new_detid_,   wrapper_uint,  clean_cols_detid) // ce_clean_detid
    .Define("ce_clean_x_unshifted",       wrapper_float, clean_cols_x)
    .Define("ce_clean_y_unshifted",       wrapper_float, clean_cols_y)
    .Define(new_z_,       wrapper_float, clean_cols_z)
    .Define(new_layer_,   wrapper_uint,  clean_cols_layer) // ce_clean_layer
    .Define(new_en_,      wrapper_float, clean_cols_en)    // ce_clean_energy

    // Define ce_clean_energy_MeV column
    .Define(new_en_MeV_, weight_energy_ce, {new_en_, new_layer_, clean_cols.back()}); //clean_cols.back() is "st" ie showertype

  if(this->showertype == SHOWERTYPE::HAD) {
    partial_process = partial_process.Define(new_ahc_en_,  clean_ahc<float>, {"ahc_hitEnergy", "ahc_hitEnergy", "ahc_hitK", clean_cols.back()})
    .Define(new_ahc_en_MeV_, weight_energy_ahc, {new_ahc_en_, clean_cols.back()});
  }


  ROOT::RDF::RNode partial_process_next = partial_process;

  /*
  The following code deals with impact and shifts.
  see the notebook https://github.com/tcuisset/HgcalClue3DClusteringPlotting/tree/main/notebooks/impact-shift-rechits.ipynb 
  for a study of all this

  Shifting positions : see https://gitlab.cern.ch/cms-hgcal-tb/TestBeam/-/issues/17
  - in data : everything is more or less misaligned (DWC, beam, detector as a whole, individual layers...)
  DWC coordinates are mirrored compared to HGCAL coordinates
  There is an Impact_Shifts.txt file (computed using muon runs) that makes the mapping between DWC coordinates and real HGCAL layer positions

  - in simulation : everything seems to be perfectly aligned
  DWC coordinates are the right way around (so different than data...)
  No need to apply Impact_Shifts.txt (everything is already centered at 0)

  The output tree has :
  - for data : you have two choices
       - either use impactX_shifted and ce_clean_x_unshifted. These two branches can be directly compared (impactX was shifted so it matches ce_clean_x)
         It does not reproduce in ce_clean_x the misalignment of the different layers, but it can directly be compared to simulation
       - either use impactX_unshifted and ce_clean_x_shifted. These 2 branches can be compared (in this case rechits positions were shifted to match extrapolated DWC tracks)
         The misalignment of layers is reproduced in hits position.
         Cannot really be compared with simulation as no simulation with misaligned layers seems to exist
  - for simulation, only one choice:
       - use impactX_unshifted and ce_clean_x_unshifted. All layers are aligned, and no shift needs to be applied
  */
  if (this->datatype == DATATYPE::DATA) {
    partial_process_next = partial_process
      /* Create new columns impactX and impactY, of type std::vector<float>, holding extrapolated DWC impact positions per layer
      (not shifted). This just packs all the myFriend.impactX_HGCal_layer_i into a vector and then multiplies it by -1
      The -1 is due to the mirrored DWC coordinates in data only (see https://gitlab.cern.ch/cms-hgcal-tb/TestBeam/-/issues/17 comment by Artur Lobanov)
      */
      .Define("impactX_unshifted", ROOT::RDF::PassAsVec<static_cast<unsigned>(detectorConstants::totalnlayers), float>([](std::vector<float> v){
        for (float& impact : v)
          impact *= -1.;
        return v;
      }), impactXcols_)
      .Define("impactY_unshifted", ROOT::RDF::PassAsVec<static_cast<unsigned>(detectorConstants::totalnlayers), float>([](std::vector<float> v){
        for (float& impact : v)
          impact *= -1.;
        return v;
      }), impactYcols_);

    /**
     * DATA ONLY (no shift in simulation)
     * Makes new column impactX_shifted : vector<float> = { -1 * impactX_HGCal_layer_i + shift[layer i] }
     * \param v branch impactX_unshifted
    */
    auto shift_impactX = [this](const std::vector<float>& v) {
      std::vector<float> newv(v.size());
      for(unsigned i=0; i<v.size(); ++i)
        newv[i] = v[i] + this->shifts_map_[i].first;
      return newv;
    };

    auto shift_impactY = [this](const std::vector<float>& v) {
      std::vector<float> newv(v.size());
      for(unsigned i=0; i<v.size(); ++i)
        newv[i] = v[i] + this->shifts_map_[i].second;
      return newv;
    };
    partial_process_next = partial_process_next
      .Define("impactX_shifted", shift_impactX, {"impactX_unshifted"})
      .Define("impactY_shifted", shift_impactY, {"impactY_unshifted"});

    savedcols_.push_back("impactX_shifted");
    savedcols_.push_back("impactY_shifted");


    /*
    DATA ONLY
    Make new column ce_clean_x_shifted which holds rechit x position shifted by shifts in Impact_Shifts.txt 
    */
    auto shift_rechitX = [this](const std::vector<float>& rechit_x, const std::vector<float>& en, const std::vector<unsigned>& l, const std::vector<unsigned>& chip, const std::vector<unsigned>& channel, const std::vector<unsigned>& module, const std::vector<float>& amplitude, const std::vector<bool>& noise_flag, const bool& st)
    {
      std::vector<float> rechit_x_cleaned = clean_ce<float>(rechit_x, en, l, chip, channel, module, amplitude, noise_flag, this->noise_map_, st);
      for(unsigned i = 0; i < rechit_x_cleaned.size(); ++i)
        rechit_x_cleaned[i] -= this->shifts_map_[l[i]].first;
      return rechit_x_cleaned;
    };
    auto shift_rechitY = [this](const std::vector<float>& rechit_y, const std::vector<float>& en, const std::vector<unsigned>& l, const std::vector<unsigned>& chip, const std::vector<unsigned>& channel, const std::vector<unsigned>& module, const std::vector<float>& amplitude, const std::vector<bool>& noise_flag, const bool& st)
    {
      std::vector<float> rechit_y_cleaned = clean_ce<float>(rechit_y, en, l, chip, channel, module, amplitude, noise_flag, this->noise_map_, st);
      for(unsigned i = 0; i < rechit_y_cleaned.size(); ++i)
        rechit_y_cleaned[i] -= this->shifts_map_[l[i]].second;
      return rechit_y_cleaned;
    };

    partial_process_next = partial_process_next
      /* Make new columns with rechits x and y positions shifted by the Shift_Impacts.txt file*/
      .Define("ce_clean_x_shifted", shift_rechitX, clean_cols_x)
      .Define("ce_clean_y_shifted", shift_rechitY, clean_cols_y);
    savedcols_.push_back("ce_clean_x_shifted");
    savedcols_.push_back("ce_clean_y_shifted");

  } else { //Monte Carlo
    /* In case of MC the impact branches are the right way around so just do nothing on the values
    */
    partial_process_next = partial_process
      .Define("impactX_unshifted", ROOT::RDF::PassAsVec<static_cast<unsigned>(detectorConstants::totalnlayers), float>([](std::vector<float> v){
        return v;
      }), impactXcols_)
      .Define("impactY_unshifted", ROOT::RDF::PassAsVec<static_cast<unsigned>(detectorConstants::totalnlayers), float>([](std::vector<float> v){
        return v;
      }), impactYcols_);
  }

  if (this->datatype == DATATYPE::MC) {
    // Save gun energy
    savedcols_.emplace_back("trueBeamEnergy");
  }
    
  if(this->datatype == DATATYPE::MC and this->showertype == SHOWERTYPE::HAD) {
    partial_process_next.Filter("ahc_energySum == 0").Snapshot(this->outdata_.tree_name.c_str(), this->outdata_.file_path.c_str(), savedcols_);
  }
  else if(this->datatype == DATATYPE::DATA or this->showertype == SHOWERTYPE::EM) {
    if(this->showertype == SHOWERTYPE::EM)
      savedcols_.erase(std::remove(savedcols_.begin(), savedcols_.end(), new_ahc_en_MeV_), savedcols_.end()); //erase-remove idiom
    partial_process_next.Snapshot(this->outdata_.tree_name.c_str(), this->outdata_.file_path.c_str(), savedcols_);
  }
}

void Selector::print_relevant_branches(const int& nrows=5, std::optional<std::string> filename)
{
  //define dataframe that owns the TTree
  ROOT::RDataFrame d(this->indata_.tree_name.c_str(), this->indata_.file_path.c_str());
  //prints TTree content according to the specified columns; a multi-threaded version is not supported
  auto display = d.Display(savedcols_, nrows);
  //write contents to file
  std::ofstream f;
  if(filename != std::nullopt)
    f.open(filename.value());
  else
    f.open("data.txt");
  f << display->AsString();
  f.close();
};

int Selector::sanity_checks(const std::string& fname) 
{
  if(fname.substr(fname.length()-5) != ".root")
    throw std::invalid_argument("The file has to be i the ROOT format.");
  return 1;
}

template std::vector<float> Selector::clean_ce(const std::vector<float>&, const std::vector<float>&, const std::vector<unsigned>&, const std::vector<unsigned>&, const std::vector<unsigned>&, const std::vector<unsigned>&, const std::vector<float>&, const std::vector<bool>&, const mapT&, const bool&);
template std::vector<unsigned> Selector::clean_ce(const std::vector<unsigned>&, const std::vector<float>&, const std::vector<unsigned>&, const std::vector<unsigned>&, const std::vector<unsigned>&, const std::vector<unsigned>&, const std::vector<float>&, const std::vector<bool>&, const mapT&, const bool&);
template std::vector<float> Selector::clean_ahc(const std::vector<float>&, const std::vector<float>&, const std::vector<int>&, const bool&);
