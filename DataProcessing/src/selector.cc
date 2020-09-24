#include "UserCode/DataProcessing/interface/selector.h"
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
  this->savedcols_ = {"event", "run", "NRechits", new_detid_, new_x_, new_y_, new_z_, new_layer_, new_en_MeV_, new_ahc_en_MeV_, "beamEnergy", new_impX_, new_impY_};
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
  std::string home( getenv("HOME") );
  std::string rel( getenv("CMSSW_VERSION") );
  std::string filename = home + "/" + rel + "/src/UserCode/DataProcessing/Noise_Map.txt";
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
  std::string home( getenv("HOME") );
  std::string rel( getenv("CMSSW_VERSION") );
  std::string filename = home + "/" + rel + "/src/UserCode/DataProcessing/Impact_Shifts.txt";
  std::ifstream in(filename);
  if(!in.is_open())
    throw std::invalid_argument("File " + filename + " could not be opened.");

  unsigned layer;
  float shiftx, shifty;
  for(unsigned i=0; i<2; ++i) //lines to skip
    in.ignore(std::numeric_limits<unsigned>::max(), '\n');
  for(unsigned i=0; i<detectorConstants::totalnlayers; ++i) {
    in >> layer >> shiftx >> shifty;
    this->shifts_map_.emplace_back( shiftx, shifty );
  }
}

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

bool Selector::common_selection(const unsigned& layer, const float& energy, const unsigned& chip, const unsigned& channel, const unsigned& module, const float& amplitude, const bool& noise_flag, const mapT& map, const bool& showertype)
{
  if(layer<1 || layer>detectorConstants::totalnlayers)
    throw std::out_of_range("Unphysical layer number: " + std::to_string(layer));
  else if(!showertype and layer>detectorConstants::nlayers_emshowers)
    return false;

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

//cleans hits in CE-E and CE-H
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

//weights the energy of the CE-E and CE-H hits
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

      if( (!st and layer>0 and layer<=detectorConstants::nlayers_emshowers)
	  or (st and layer>0 and layer<=detectorConstants::totalnlayers) )
	{
	  if(st and layer>detectorConstants::nlayers_emshowers)
	    weight = detectorConstants::globalWeightCEH;
	  else
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

  if(!st)
    weight = 0.f;

  for(unsigned i=0; i<nhits; ++i)
    en_weighted.push_back( weight * en[i] );
      
  return en_weighted;  
}

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

  f_had = new TFile(this->indata_.file_path.c_str());
  t_had1 = static_cast<TTree*>( f_had->Get(this->indata_.tree_name.c_str()) ); 
  t_had2 = static_cast<TTree*>( f_had->Get(this->indata_.tree_name_friend.c_str()) ); 
  t_had1->AddFriend(t_had2, "myFriend");
  d = new ROOT::RDataFrame(*t_had1);

  ROOT::EnableImplicitMT( ncpus_ ); //enable parallelism
  ROOT::Detail::RDF::ColumnNames_t clean_cols = {"rechit_energy", "rechit_layer", "rechit_chip", "rechit_channel", "rechit_module", "rechit_amplitudeHigh", "rechit_noise_flag", "st"};
  ROOT::Detail::RDF::ColumnNames_t clean_cols_detid = clean_cols; clean_cols_detid.insert( clean_cols_detid.begin(), "rechit_detid");
  ROOT::Detail::RDF::ColumnNames_t clean_cols_x = clean_cols; clean_cols_x.insert( clean_cols_x.begin(), "rechit_x");
  ROOT::Detail::RDF::ColumnNames_t clean_cols_y = clean_cols; clean_cols_y.insert( clean_cols_y.begin(), "rechit_y");
  ROOT::Detail::RDF::ColumnNames_t clean_cols_z = clean_cols; clean_cols_z.insert( clean_cols_z.begin(), "rechit_z");
  ROOT::Detail::RDF::ColumnNames_t clean_cols_layer = clean_cols; clean_cols_layer.insert( clean_cols_layer.begin(), "rechit_layer");
  ROOT::Detail::RDF::ColumnNames_t clean_cols_en = clean_cols; clean_cols_en.insert( clean_cols_en.begin(), "rechit_energy");
  
  //wrappers required to circumvent the static-ness of the wrapped methods (required by RDataFrame) while passing a class variable (this->noise_map_)
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

  auto shift_impactX = [this](const std::vector<float>& v) {
    std::vector<float> newv(v.size());
    for(unsigned i=0; i<v.size(); ++i)
      newv[i] = (-1.f * v[i]) + this->shifts_map_[i].first;
    return newv;
  };

  auto shift_impactY = [this](const std::vector<float>& v) {
    std::vector<float> newv(v.size());
    for(unsigned i=0; i<v.size(); ++i)
      newv[i] = (-1.f * v[i]) + this->shifts_map_[i].second;
    return newv;
  };

  std::string filter_str = "true", define_str = "return false;";
  if(this->showertype == SHOWERTYPE::HAD) {
    filter_str = "myFriend.dwcReferenceType>=13 && myFriend.trackChi2_X<10 && myFriend.trackChi2_Y<10";
    define_str = "return true;";
  }

  auto partial_process = d->Filter(filter_str.c_str())
    .Filter(ROOT::RDF::PassAsVec<static_cast<unsigned>(2*detectorConstants::totalnlayers), float>(remove_missing_dwc),
	    impactcols_)
    .Define(new_impX_, ROOT::RDF::PassAsVec<static_cast<unsigned>(detectorConstants::totalnlayers), float>(shift_impactX), impactXcols_)
    .Define(new_impY_, ROOT::RDF::PassAsVec<static_cast<unsigned>(detectorConstants::totalnlayers), float>(shift_impactY), impactYcols_)
    .Define(clean_cols.back(), define_str) //showertype: em or had
    .Define(new_detid_,   wrapper_uint,  clean_cols_detid)
    .Define(new_x_,       wrapper_float, clean_cols_x)
    .Define(new_y_,       wrapper_float, clean_cols_y)
    .Define(new_z_,       wrapper_float, clean_cols_z)
    .Define(new_layer_,   wrapper_uint,  clean_cols_layer)
    .Define(new_en_,      wrapper_float, clean_cols_en)
    .Define(new_en_MeV_, weight_energy_ce, {new_en_, new_layer_, clean_cols.back()});

  if(this->showertype == SHOWERTYPE::HAD) {
    partial_process = partial_process.Define(new_ahc_en_,  clean_ahc<float>, {"ahc_hitEnergy", "ahc_hitEnergy", "ahc_hitK", clean_cols.back()})
    .Define(new_ahc_en_MeV_, weight_energy_ahc, {new_ahc_en_, clean_cols.back()});
  }
    
  if(this->datatype == DATATYPE::MC and this->showertype == SHOWERTYPE::HAD) {
    partial_process.Filter("ahc_energySum == 0").Snapshot(this->outdata_.tree_name.c_str(), this->outdata_.file_path.c_str(), savedcols_);
  }
  else if(this->datatype == DATATYPE::DATA or this->showertype == SHOWERTYPE::EM) {
    if(this->showertype == SHOWERTYPE::EM)
      savedcols_.erase(std::remove(savedcols_.begin(), savedcols_.end(), new_ahc_en_MeV_), savedcols_.end()); //erase-remove idiom
    partial_process.Snapshot(this->outdata_.tree_name.c_str(), this->outdata_.file_path.c_str(), savedcols_);
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
