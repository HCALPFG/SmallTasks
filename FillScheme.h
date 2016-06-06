//Reads in a Bunch Fill file and returns the scehme as a vector<bool>
vector<bool> readFillScheme(TString indir, int fill){

  vector<bool> fillscheme;
  ifstream file(indir+"BunchFill"+to_string(fill)+".txt");
  //Check that the file exists
  if(!file.is_open()){ cout<<"\n\e[31m[HF Timing]\e[0m ERROR: Fill scheme file does not exist for fill "<<fill<<"!"<<endl; exit(0);}
  string str;
  while(getline(file,str))
    fillscheme.push_back(stoi(str));

  return fillscheme;
}

//Takes in the fillmap, run list, and loops over and reads Bunch Fill files from indir
vector<vector<bool> > getFillSchemes(map<unsigned int,int> fillmap, vector<unsigned int> runs, TString indir){
  cout<<"[HF Timing] Getting fill schemes...\n"<<endl;

  vector<vector<bool> > fillSchemes;

  int tmp_fill = -1;
  vector<bool> tmp_scheme;
  for(unsigned int irun=0; irun<runs.size(); irun++){
    int curr_fill = fillmap[runs[irun]];
    //Don't need to read in file again for run in same fill
    if(curr_fill!=tmp_fill){
      tmp_scheme = readFillScheme(indir, curr_fill);
      fillSchemes.push_back(tmp_scheme);
      tmp_fill = curr_fill;
    }
    else{
      fillSchemes.push_back(tmp_scheme);
    }
  }
  return fillSchemes;
}

//Takes run-list vector and vector of fillschemes and finds the bxs according to bxtype
map<unsigned int, vector<int> > selectBXs(vector<unsigned int> runs, vector<vector<bool> > schemes, TString bxtype){
  cout<<"[HF Timing] Selecting "<<bxtype<<" bxs..."<<endl;
  
  if(bxtype!="iso" && bxtype!="first" && bxtype!="last" && bxtype!="middle" && bxtype!="noniso"){
    cout<<"\n\e[31m[HF Timing]\e[0m WARNING: Not valid bx type!"<<endl;
    exit(0);
  }
  
  // Map betwen run and vector of selected bxs
  map<unsigned int, vector<int> > selectedbx_map;
  for(unsigned int irun=0; irun<runs.size(); irun++){
    vector<int> tmp_bxs;
    
    for(unsigned int ibx=1; ibx<3564; ibx++){ //Start at 1
      if(schemes[irun][ibx]){ //veto empty bxs
	
	// Not most efficient way, but clearer code
	if(bxtype=="iso" && !schemes[irun][ibx-1] && !schemes[irun][ibx+1])
	  tmp_bxs.push_back(ibx);
	else if(bxtype=="first" && !schemes[irun][ibx-1] && schemes[irun][ibx+1])
	  tmp_bxs.push_back(ibx);
	else if(bxtype=="last" && schemes[irun][ibx-1] && !schemes[irun][ibx+1])
	  tmp_bxs.push_back(ibx);
	else if(bxtype=="middle" && schemes[irun][ibx-1] && schemes[irun][ibx+1])
	  tmp_bxs.push_back(ibx);
	else if(bxtype=="noniso" && !(!schemes[irun][ibx-1] && !schemes[irun][ibx+1]))
	  tmp_bxs.push_back(ibx);
      }
    }
    if(tmp_bxs.size()==0) cout<<"\e[31m[HF Timing]\e[0m WARNING: Run "<<runs[irun]<<" does not have any "<<bxtype<<" bunches!"<<endl;
    selectedbx_map[runs[irun]] = tmp_bxs;
  }
  return selectedbx_map;
}

//Takes in the fillmap and selected bxs and prints them to outdir in a txt file. Can be toggled on and off with print bool.
void printSelectedBXs(bool print, map<unsigned int,int> fillmap, map<unsigned int, vector<int> > selectedbxs, TString outdir){

  if(print){
    ofstream file(outdir+"selectedbxs.txt");
    for(map<unsigned int,vector<int> >::iterator imap=selectedbxs.begin(); imap!=selectedbxs.end(); ++imap){
      file<<"Fill: "<<fillmap[(*imap).first]<<", Run: "<<(*imap).first<<", Nbxs: "<<(*imap).second.size()<<endl;
      for(unsigned int i=0;i<(*imap).second.size();i++){
	file<<"bx = "<<(*imap).second[i]<<endl;
      }
    }
    file.close();
  }
}

//Pass vector of selected bxs and the bx to check
bool isSelectedBX(vector<int> selectedBXs, int bx){

  bool isSelected = false;
  if(find(selectedBXs.begin(), selectedBXs.end(), bx) != selectedBXs.end())
    isSelected = true;

  return isSelected;
}
