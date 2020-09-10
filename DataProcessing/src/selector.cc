#include "UserCode/DataProcessing/interface/selector.h"

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
}

Selector::~Selector()
{
}

bool Selector::common_selection(const unsigned& layer, const float& energy, const unsigned& chip, const unsigned& channel, const unsigned& module, const bool& noise_flag)
{
  if(layer<1 || layer>detectorConstants::totalnlayers)
    throw std::out_of_range("Unphysical layer number: " + std::to_string(layer));
  else if(this->showertype == SHOWERTYPE::EM and layer>detectorConstants::nlayers_emshowers && layer<=detectorConstants::totalnlayers)
    return false;
  else if(this->showertype == SHOWERTYPE::HAD and layer>=36 and layer<=37)
    return false;

  if(noise_flag==1)
    return false;

  if(energy<=0.5) //MIP noise cut
    return false;
  if( (chip==0 and layer==1) or (chip==0 and channel==44) or (chip==3 and channel==22) or (chip==3 and channel==28) ) //mask noisy chips
    return false;
  if( chip==1 and layer==36 and channel==16 and module==39 )
    return false;
  return true;
}

void Selector::select_relevant_branches()
{  
  auto clean_float_arrays = [&](std::vector<float> var, std::vector<float> en, std::vector<unsigned> l, std::vector<unsigned> chip, std::vector<unsigned> channel, std::vector<unsigned> module, std::vector<bool> noise_flag) {
      unsigned layer = 0;
      float energy = 0.f;
      size_t ensize = en.size();
      std::vector<float> var_clean;
      var_clean.reserve(ensize);

      for(unsigned i=0; i<ensize; ++i) {
	energy = en[i];
	layer = l[i];
	if( common_selection(layer, energy, chip[i], channel[i], module[i], noise_flag[i]) )
	  var_clean.push_back( var[i] );
      }
      return var_clean;  
    };

  auto clean_detids = [&](std::vector<unsigned> detid, std::vector<float> en, std::vector<unsigned> l, std::vector<unsigned> chip, std::vector<unsigned> channel, std::vector<unsigned> module, std::vector<bool> noise_flag) {
      unsigned layer = 0;
      float energy = 0.f;
      size_t ensize = en.size();
      std::vector<unsigned> detid_clean;
      detid_clean.reserve(ensize);

      for(unsigned i=0; i<ensize; ++i)
	{
	  energy = en[i];
	  layer = l[i];
	  if( common_selection(layer, energy, chip[i], channel[i], module[i], noise_flag[i]) )
	    detid_clean.push_back( detid[i] );
	}
      return detid_clean;  
    };

  auto clean_layers = [&](std::vector<float> en, std::vector<unsigned> l, std::vector<unsigned> chip, std::vector<unsigned> channel, std::vector<unsigned> module, std::vector<bool> noise_flag) {
      unsigned layer = 0;
      float energy = 0.f;
      size_t ensize = en.size();
      std::vector<unsigned> layer_clean;
      layer_clean.reserve(ensize);

      for(unsigned i=0; i<ensize; ++i) {
	energy = en[i];
	layer = l[i];
	if( common_selection(layer, energy, chip[i], channel[i], module[i], noise_flag[i]) )
	  layer_clean.push_back( layer );
      }
      return layer_clean;  
    };
      
  auto weight_and_clean_energy = [&](std::vector<float> en, std::vector<unsigned> l, std::vector<unsigned> chip, std::vector<unsigned> channel, std::vector<unsigned> module, std::vector<bool> noise_flag) {
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

	  if( (this->showertype==SHOWERTYPE::EM and layer>0 and layer<=detectorConstants::nlayers_emshowers)
	      or (this->showertype==SHOWERTYPE::HAD and layer>0 and layer<=detectorConstants::totalnlayers) )
	    {
	      if(this->showertype==SHOWERTYPE::HAD and layer>detectorConstants::nlayers_emshowers) //remove as soon as the extra dEdX are known
		weight = 1.;
	      else
		{
		  //remove this bit once the FH weights are established
		  float XXXXweight;
		  if(layer > detectorConstants::nlayers_emshowers)
		    XXXXweight = 1.f;
		  else
		    XXXXweight = detectorConstants::dEdX.at(layer-1);
		  ///////////////////////////////////////////////////

		  weight = XXXXweight;
		}
	    }
	  else
	    continue;

	  if( common_selection(layer, energy, chip[i], channel[i], module[i], noise_flag[i]) )
	    en_weighted.push_back( energy * weight );
	}
      return en_weighted;  
    };

  //noise rejection cuts
  auto noise_rejection = [&](const std::vector<unsigned> l) {
    unsigned counter_cee=0, counter_ceh=0;

    for(unsigned i=0; i<l.size(); ++i)
      {
	if(l[i]>0 and l[i]<detectorConstants::nlayers_emshowers)
	  ++counter_cee;
	if(this->showertype==SHOWERTYPE::HAD) {
	  if( (l[i]>detectorConstants::nlayers_emshowers and l[i]<36) or (l[i]>37 and l[i]<=detectorConstants::totalnlayers) )
	    ++counter_ceh;
	}
      }

    if( counter_cee>18
	or (this->showertype==SHOWERTYPE::HAD and counter_ceh>24) ) //6 * #sigmas, where 6 is the approx. average of noise in the chips
      return true;
    return false;
  };

    //missing s1/s25 cut
    /*
      auto muon_veto = [&](const std::vector<float>& en, const std::vector<unsigned>& l) {
      if(this->showertype == SHOWERTYPE::EM)
      return true;
      else if(this->showertype == SHOWERTYPE::HAD)
      {
      float recoen_cee=0.f, recoen_ceh=0.f;
      for(unsigned i=0; i<en.size(); ++i) 
      {
      if(l[i]<=detectorConstants::nlayers_emshowers)
      recoen_cee += en[i];
      else
      recoen_ceh += en[i];
      }
      if(recoen_cee<100 and recoen_ceh<60)
      return false;
      }
      return true;
      };
    */
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

    ROOT::EnableImplicitMT( ncpus_ ); //enable parallelism

    std::string filter_str = "true";
    if(this->showertype == SHOWERTYPE::HAD)
      filter_str = "myFriend.dwcReferenceType && myFriend.trackChi2_X<10 && myFriend.trackChi2_Y<10";
	
    auto partial_process = d->Filter(noise_rejection, {"rechit_layer"})
    //.Filter(muon_veto, {"rechit_energy"})
    .Filter(filter_str.c_str())
    .Define(this->newcol_clean_detid_,     clean_detids,              {"rechit_detid", "rechit_energy", "rechit_layer", "rechit_chip", "rechit_channel", "rechit_module", "rechit_noise_flag"})
    .Define(this->newcol_clean_x_,         clean_float_arrays,        {"rechit_x",     "rechit_energy", "rechit_layer", "rechit_chip", "rechit_channel", "rechit_module", "rechit_noise_flag" })
    .Define(this->newcol_clean_y_,         clean_float_arrays,        {"rechit_y",     "rechit_energy", "rechit_layer", "rechit_chip", "rechit_channel", "rechit_module", "rechit_noise_flag"})
    .Define(this->newcol_clean_z_,         clean_float_arrays,        {"rechit_z",     "rechit_energy", "rechit_layer", "rechit_chip", "rechit_channel", "rechit_module", "rechit_noise_flag"})
    .Define(this->newcol_clean_layer_,     clean_layers,              {"rechit_energy", "rechit_layer", "rechit_chip", "rechit_channel", "rechit_module", "rechit_noise_flag"})
    .Define(this->newcol_clean_energy_,    weight_and_clean_energy,   {"rechit_energy", "rechit_layer", "rechit_chip", "rechit_channel", "rechit_module", "rechit_noise_flag"});

    if(this->datatype == DATATYPE::MC) {
      std::cout << "snap mc" << std::endl;
      partial_process.Filter("ahc_energySum == 0").Snapshot(this->outdata_.tree_name.c_str(), this->outdata_.file_path.c_str(), this->cols_);
    }
    else if(this->datatype == DATATYPE::DATA) {
      std::cout << "snap" << std::endl;
      partial_process.Snapshot(this->outdata_.tree_name.c_str(), this->outdata_.file_path.c_str(), this->cols_);
    }
};

void Selector::print_relevant_branches(const int& nrows=5, std::optional<std::string> filename)
{
  //define dataframe that owns the TTree
  ROOT::RDataFrame d(this->indata_.tree_name.c_str(), this->indata_.file_path.c_str());
  //prints TTree content according to the specified columns; a multi-threaded version is not supported
  auto display = d.Display(cols_, nrows);
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
