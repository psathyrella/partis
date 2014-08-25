//----------------------------------------------------------------------------------------
void set_bins(vector<double> values, int n_bins, bool is_log_x, TString var_type, double *xbins, map<TString,int> string_map) {
  if (is_log_x) {
    double log_xmin = log(powf(10.0f, floorf(log10f(values[0]))));  // round down to the nearest power of 10
    double log_xmax = log(powf(10.0f, ceilf(log10f(values[values.size()-1])))); // round up to the nearest power of 10
    double log_dx = (log_xmax - log_xmin) / n_bins;
    log_xmin -= 0.1*log_dx;  // expand a little to avoid overflows
    log_xmax += 0.1*log_dx;
    log_dx = (log_xmax - log_xmin) / n_bins;
    for (int ib=0; ib<=n_bins; ib++) {
      double low_edge = exp(log_xmin + ib * log_dx);
      xbins[ib] = low_edge;
    }
  } else {
    double dx = (values[values.size()-1] - values[0]) / n_bins;
    double xmin = values[0] - 0.1*dx;  // expand a little to avoid overflows
    double xmax = values[values.size()-1] + 0.1*dx;  // expand a little to avoid overflows
    if (var_type=="string") {
      xmin = 0.5;
      xmax = string_map.size() + 0.5;
    }
    dx = (xmax - xmin) / n_bins;  // then recalculate dx
    for (int ib=0; ib<=n_bins; ib++) {
      xbins[ib] = xmin + ib * dx;
    }
  }
}
//----------------------------------------------------------------------------------------
TH1F make_hist(TString infname, TString data_type, TString var_type, TString log, map<TString,int> string_map,
	       int n_bins=30, float xmin_force=0.0, float xmax_force=0.0, int n_dummy_columns=0, int category_column=-1) {

  if (var_type=="string")
    n_bins = string_map.size();
  double xbins[n_bins+1];  // NOTE has to be n_bins + 1

  ifstream ifs(infname);
  if(!ifs.is_open()) {
    cout << "    " << infname << " d.n.e." << endl;
    return TH1F();
  }
  string line;
  vector<double> values;
  while(getline(ifs,line)) {
    stringstream ss(line);

    // skip some columns
    int idummy(0);
    string sdummy;
    while (idummy < n_dummy_columns) {
      ss >> sdummy;
      ++idummy;
    }

    if (category_column >= 0) {  // require that the column immediately before the data column has value <category_column>
      int category;
      ss >> category;
      if (category != category_column)
	continue;
    }

    if (var_type=="double") {
      float value;
      ss >> value;
      values.push_back(value);
    } else if (var_type=="string") {
      string value;
      ss >> value;
      if (string_map.find(value) == string_map.end()) {
	cout << value << " not found!" << endl;
	assert(0);
      }
      values.push_back(double(string_map[value] + 1));  // string_map is a crappy name, but I can't come up with anything better
                                                        // NOTE add one to get to 1-base indexing for hist bins
   } else {
      assert(0);
    }
  }
  if (values.size() == 0)
    return TH1F();
  sort(values.begin(), values.end());

  TH1F hist;
  if (xmin_force == xmax_force) {  // if boundaries aren't set explicitly, work them out dynamically
    set_bins(values, n_bins, log.Contains("x"), var_type, xbins, string_map);
    hist = TH1F("h"+data_type, "", n_bins, xbins);
  } else {
    hist = TH1F("h"+data_type, "", n_bins, xmin_force, xmax_force);
  }
  hist.Sumw2();
  for (unsigned iv=0; iv<values.size(); iv++) {
    hist.Fill(values[iv]);
  }

  cout << "  " << values.size() << " (mean " << hist.GetMean() << ") values in " << infname << endl;

  // make sure there's no overflows
  if(hist.GetBinContent(0) != 0 || hist.GetBinContent(hist.GetNbinsX()+1) != 0) {
    cout << infname << endl;
    for (unsigned iv=0; iv<values.size(); iv++) {
      cout
	<< setw(12) << values[iv];
    }
    cout << endl;
    for (int ib=0; ib<hist.GetNbinsX()+2; ib++) {
      cout
      	<< setw(12) << ib
      	<< setw(12) << hist.GetBinLowEdge(ib)
      	<< setw(12) << hist.GetBinContent(ib)
      	<< endl;
    }
    // assert(0);
    cout << "WARNING you got overflows!" << endl;
  }

  hist.Scale(1./hist.Integral());
  
  return hist;
}
