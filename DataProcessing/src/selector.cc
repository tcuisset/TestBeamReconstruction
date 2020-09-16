#include "UserCode/DataProcessing/interface/selector.h"

Selector::Selector(const std::string& in_file_path, const std::string& out_file_path,
		   const std::string& datatype_, const std::string& showertype_, const int& beam_energy_,
		   std::optional<std::string> in_tree_name, std::optional<std::string> in_tree_name_friend, std::optional<std::string> out_tree_name)
{
  std::cout << "check1" << std::endl;
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
  std::cout << "check2" << std::endl;
  this->load_noise_values();
  std::cout << "check3" << std::endl;
}

Selector::~Selector()
{
}

void Selector::load_noise_values()
{
  std::cout <<  "load" << std::endl;
  std::string home( getenv("HOME") );
  std::cout <<  "home" << std::endl;
  std::string rel( getenv("CMSSW_VERSION") );
  std::cout <<  "version" << std::endl;
  std::string filename = home + "/" + rel + "/src/UserCode/DataProcessing/Noise_Map.txt";
  std::cout << "Filename for noise map: " << filename  << std::endl;
  std::ifstream in(filename);
  if(!in.is_open())
    throw std::invalid_argument("File " + filename + " could not be opened.");

  std::pair<int,int> mod_chip;
  std::pair<std::pair<int,int>, float> temp;
  unsigned layer, mod_id, mod_pos, chip;
  float noise;
  while(in >> layer >> mod_id >> mod_pos >> chip >> noise)
    {
      mod_chip = std::make_pair(mod_id, chip);
      temp = std::make_pair(mod_chip, noise);
      this->noise_map_.insert(temp);
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
  else if(!showertype and layer>detectorConstants::nlayers_emshowers && layer<=detectorConstants::totalnlayers)
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

template<typename T>
std::vector<T> Selector::clean_rechits(const std::vector<T>& var, const std::vector<float>& en, const std::vector<unsigned>& l, const std::vector<unsigned>& chip, const std::vector<unsigned>& channel, const std::vector<unsigned>& module, const std::vector<float>& amplitude, const std::vector<bool>& noise_flag, const mapT& map, const bool& st)
{
  unsigned layer = 0;
  float energy = 0.f;
  size_t ensize = en.size();
  std::vector<T> var_clean;
  var_clean.reserve(ensize);

  for(unsigned i=0; i<ensize; ++i) {
    energy = en[i];
    layer = l[i];
    if( common_selection(layer, energy, chip[i], channel[i], module[i], amplitude[i], noise_flag[i], map, st) )
      var_clean.push_back( var[i] );
  }
  return var_clean;  
}

template<typename T>
std::vector<T> Selector::clean_ahc(const std::vector<T>& var, const std::vector<float>& en, const std::vector<int>& l, const bool& st)
{
  unsigned layer = 0;
  size_t size = var.size();
  std::vector<T> var_clean;
  var_clean.reserve(size);

  if(st)
    {
      for(unsigned i=0; i<size; ++i) {
	layer = l[i];
	if(layer != 38) //mask AHCAL layer #38
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

std::vector<float> Selector::weight_energy(const std::vector<float>& en, const std::vector<unsigned>& l, const float& ahc_en, const bool& st) 
{
  unsigned layer = 0;
  float energy = 0.f;
  float weight = 1.f;
  size_t ensize = en.size();
  std::vector<float> en_weighted;
  en_weighted.reserve(ensize); //maximum possible size

  for(unsigned i=0; i<ensize; ++i)
    {
      energy = en[i];
      layer = l[i];

      if( (!st and layer>0 and layer<=detectorConstants::nlayers_emshowers)
	  or (st and layer>0 and layer<=detectorConstants::totalnlayers) )
	{
	  if(st and layer>detectorConstants::nlayers_emshowers) //remove as soon as the extra dEdX are known
	    weight = detectorConstants::globalWeightCEH;
	  else
	    weight = detectorConstants::dEdX.at(layer-1);
	}
      else
	continue;

      if(layer>detectorConstants::nlayers_emshowers)
	en_weighted.push_back( weight * (energy + 0.4*ahc_en) );
      else
	en_weighted.push_back( energy * weight );
 
    }
  return en_weighted;  
}

/*
unsigned calculate_shower_start(const bool& st)
{
  if(!st) //positron showers do not required this calculation
    return 0;

  for(unsigned il=0; il<detectorConstants::totalnlayers; ++il)
    {
      
    }
  return 0;
}
*/

void Selector::select_relevant_branches()
{
  std::cout << "check4" << std::endl;
  ROOT::RDataFrame *d;
  TFile *f_had = nullptr;
  TTree *t_had1 = nullptr;
  TTree *t_had2 = nullptr;
  if(showertype == SHOWERTYPE::HAD)
    {
      f_had = new TFile(this->indata_.file_path.c_str());
      t_had1 = static_cast<TTree*>( f_had->Get(this->indata_.tree_name.c_str()) ); 
      t_had2 = static_cast<TTree*>( f_had->Get(this->indata_.tree_name_friend.c_str()) ); 
      t_had1->AddFriend(t_had2, "myFriend");
      d = new ROOT::RDataFrame(*t_had1);
    }
  else if(showertype == SHOWERTYPE::EM)
    d = new ROOT::RDataFrame(this->indata_.tree_name.c_str(), this->indata_.file_path.c_str());

  std::cout << "check5" << std::endl;
  ROOT::EnableImplicitMT( ncpus_ ); //enable parallelism
  ROOT::Detail::RDF::ColumnNames_t clean_cols = {"rechit_energy", "rechit_layer", "rechit_chip", "rechit_channel", "rechit_module", "rechit_amplitudeHigh", "rechit_noise_flag", "st"};
  ROOT::Detail::RDF::ColumnNames_t clean_cols_detid = clean_cols; clean_cols_detid.insert( clean_cols_detid.begin(), "rechit_detid");
  ROOT::Detail::RDF::ColumnNames_t clean_cols_x = clean_cols; clean_cols_x.insert( clean_cols_x.begin(), "rechit_x");
  ROOT::Detail::RDF::ColumnNames_t clean_cols_y = clean_cols; clean_cols_y.insert( clean_cols_y.begin(), "rechit_y");
  ROOT::Detail::RDF::ColumnNames_t clean_cols_z = clean_cols; clean_cols_z.insert( clean_cols_z.begin(), "rechit_z");
  ROOT::Detail::RDF::ColumnNames_t clean_cols_layer = clean_cols; clean_cols_layer.insert( clean_cols_layer.begin(), "rechit_layer");
  ROOT::Detail::RDF::ColumnNames_t clean_cols_en = clean_cols; clean_cols_en.insert( clean_cols_en.begin(), "rechit_energy");
  
  //wrappers required to circumvent the static-ness of the wrapped methods (required by RDataFrame) while passing a class variable (this->noise_map_)
  auto wrapper_float = [&](const std::vector<float>& var, const std::vector<float>& en, const std::vector<unsigned>& l, const std::vector<unsigned>& chip, const std::vector<unsigned>& channel, const std::vector<unsigned>& module, const std::vector<float>& amplitude, const std::vector<bool>& noise_flag, const bool& st)
			{
			  std::vector<float> v = clean_rechits<float>(var, en, l, chip, channel, module, amplitude, noise_flag, this->noise_map_, st);
			  return v;
			};
  auto wrapper_uint = [&](const std::vector<unsigned>& var, const std::vector<float>& en, const std::vector<unsigned>& l, const std::vector<unsigned>& chip, const std::vector<unsigned>& channel, const std::vector<unsigned>& module, const std::vector<float>& amplitude, const std::vector<bool>& noise_flag, const bool& st)
			{
			  std::vector<unsigned> v = clean_rechits<unsigned>(var, en, l, chip, channel, module, amplitude, noise_flag, this->noise_map_, st);
			  return v;
			};

  std::string filter_str = "true", define_str = "return false;";
  if(this->showertype == SHOWERTYPE::HAD) {
    filter_str = "myFriend.dwcReferenceType && myFriend.trackChi2_X<10 && myFriend.trackChi2_Y<10";
    define_str == "return true;";
  }
  std::cout << "check6" << std::endl;
  auto partial_process = d->Filter(filter_str.c_str())
    .Define(clean_cols.back(), define_str) //showertype: em or had
    //.Define("shower_start", calculate_shower_start, {})
    .Define(new_detid_,   wrapper_uint,  clean_cols_detid)
    .Define(new_x_,       wrapper_float, clean_cols_x)
    .Define(new_y_,       wrapper_float, clean_cols_y)
    .Define(new_z_,       wrapper_float, clean_cols_z)
    .Define(new_layer_,   wrapper_uint,  clean_cols_layer)
    .Define(new_en_,      wrapper_float, clean_cols_en)
    .Define(new_ahc_en_,  clean_ahc<float>, {"ahc_hitEnergy", "ahc_hitEnergy", "ahc_hitK", clean_cols.back()})
    .Define(new_ahc_en_sum_, ahc_energy_sum, {new_ahc_en_})
    .Define(new_en_MeV_, weight_energy, {new_en_, new_layer_, new_ahc_en_sum_, clean_cols.back()});
    
  if(this->datatype == DATATYPE::MC) {
    partial_process.Filter("ahc_energySum == 0").Snapshot(this->outdata_.tree_name.c_str(), this->outdata_.file_path.c_str(), savedcols_);
  }
  else if(this->datatype == DATATYPE::DATA) {
    partial_process.Snapshot(this->outdata_.tree_name.c_str(), this->outdata_.file_path.c_str(), savedcols_);
  }
  std::cout << "check7" << std::endl;
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

template std::vector<float> Selector::clean_rechits(const std::vector<float>&, const std::vector<float>&, const std::vector<unsigned>&, const std::vector<unsigned>&, const std::vector<unsigned>&, const std::vector<unsigned>&, const std::vector<float>&, const std::vector<bool>&, const mapT&, const bool&);
template std::vector<unsigned> Selector::clean_rechits(const std::vector<unsigned>&, const std::vector<float>&, const std::vector<unsigned>&, const std::vector<unsigned>&, const std::vector<unsigned>&, const std::vector<unsigned>&, const std::vector<float>&, const std::vector<bool>&, const mapT&, const bool&);
template std::vector<float> Selector::clean_ahc(const std::vector<float>&, const std::vector<float>&, const std::vector<int>&, const bool&);
