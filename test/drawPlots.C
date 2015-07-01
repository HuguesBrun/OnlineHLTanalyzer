TFile *fileDATA = new TFile("theHistosFile_Data.root");
TFile *fileMC = new TFile("theHistosFile_MC.root");

float findMaximum(float a, float b){
    if ( a>b) return a;
    else return b;
}

void draw1Dplot(TString theName, TString Xaxis, TString Yaxis){
    TH1F *histoDATA = (TH1F*) fileDATA->Get(theName);
    TH1F *histoMC = (TH1F*) fileMC->Get(theName);
    TH1F *histoMC2 = (TH1F*) fileMC->Get(theName);
    histoMC2->Sumw2();
    
    float scale = 1.0*histoDATA->GetEntries()/histoMC->GetEntries();
    cout << "scale=" << scale << endl;
    histoMC->Scale(scale);
    histoMC2->Scale(scale);
    
    TCanvas *c0 = new TCanvas("coucou","c0",600,600);
    c0->SetFillColor(0);

    float theMaximum = 5./4*findMaximum(histoDATA->GetMaximum(), histoMC->GetMaximum());
    
    histoMC->GetXaxis()->SetTitle(Xaxis);
    histoMC->GetYaxis()->SetTitle(Yaxis);
    histoMC->SetFillColor(kGreen-7);
    histoMC2->SetFillColor(kGray+3);
    histoMC2->SetFillStyle(3001);
    histoMC->SetMaximum(theMaximum);
    histoMC->SetMinimum(0);
    histoMC->Draw();
    histoMC2->Draw("same:E2");
    
    histoDATA->Draw("E1:same");
    c0->Print("plots/"+theName+".png");
    
    
}

void draw1DdataPlot(TString theName, TString Xaxis, TString Yaxis){
    TH1F *histoDATA = (TH1F*) fileDATA->Get(theName);
    
    
    TCanvas *c0 = new TCanvas("coucou","c0",600,600);
    c0->SetFillColor(0);
    
    histoDATA->SetMinimum(0);
    histoDATA->GetXaxis()->SetTitle(Xaxis);
    histoDATA->GetYaxis()->SetTitle(Yaxis);

    
    histoDATA->Draw("E1:same");
    c0->Print("plots/"+theName+".png");
    
    
}

void draw2DdataPlot(TString theName, TString Xaxis, TString Yaxis){
    TH1F *histoDATA = (TH1F*) fileDATA->Get(theName);
    
    
    TCanvas *c0 = new TCanvas("coucou","c0",1000,600);
    c0->SetFillColor(0);
    
    histoDATA->SetMinimum(0);
    histoDATA->GetXaxis()->SetTitle(Xaxis);
    histoDATA->GetYaxis()->SetTitle(Yaxis);
    
    
    histoDATA->Draw("colz");
    c0->Print("plots/2D"+theName+".png");

}
    
void drawPlots(){
    gStyle->SetOptStat(0);
    draw1Dplot("ptL2","p_{T} L2","#");
    draw1Dplot("etaL2","#eta L2","#");
    draw1Dplot("phiL2","#phi L2","#");
    draw1Dplot("dxyL2","dxy L2","#");
    draw1Dplot("dzL2","dz L2","#");
    draw1Dplot("nbStationsWithHits","nbStationsWithHits L2","#");
    draw1Dplot("nbStationsWithHitsCSC","nbStationsWithHitsCSC L2","#");
    draw1Dplot("nbStationsWithHitsDT","nbStationsWithHitsDT L2","#");
    draw1Dplot("nbValidHits","nbValidHits L2","#");
    
    
    draw1Dplot("ptL3","p_{T} L3","#");
    draw1Dplot("etaL3","#eta L3","#");
    draw1Dplot("phiL3","#phi L3","#");
    draw1Dplot("dxyL3","dxy L3","#");
    draw1Dplot("dzL3","dz L3","#");
    draw1Dplot("drL3","dr L3","#");
    draw1Dplot("Chi2L3","Chi2 L3","#");
    
    
    draw1DdataPlot("ptL2","p_{T} L2","#");
    draw1DdataPlot("etaL2highBins","#eta L2","#");
    draw1DdataPlot("phiL2highBin","#phi L2","#");
    draw1DdataPlot("dxyL2highBin","dxy L2","#");
    draw1DdataPlot("dzL2highBin","dz L2","#");

    
    
    draw1DdataPlot("ptL3","p_{T} L3","#");
    draw1DdataPlot("etaL3highBins","#eta L3","#");
    draw1DdataPlot("phiL3highBins","#phi L3","#");
    draw1DdataPlot("dxyL3highBins","dxy L3","#");
    draw1DdataPlot("dzL3highBins","dz L3","#");
    draw1DdataPlot("drL3highBins","dr L3","#");
    draw1DdataPlot("Chi2L3highBins","Chi2 L3","#");
    
    draw2DdataPlot("L2etaPhi","L2 #phi","L2 #eta");
    draw2DdataPlot("L3etaPhi","L3 #phi","L3 #eta");

    
}