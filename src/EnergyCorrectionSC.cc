
/**
 *  @file   MarlinPandoraSC/src/EnergyCorrectionSC.cc
 * 
 *  @brief  Implementation of the SC energy correction function algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "EnergyCorrectionSC.h"
#include "PandoraPFANewProcessor.h"

//#include "TFile.h"
//#include "TH1F.h"
//#include "TTree.h"

using namespace pandora;

EnergyCorrectionSC::EnergyCorrectionSC() :
  m_SCEnergyConstants1(NULL),
  m_SCEnergyConstants2(NULL)
{
  //fCluster = new TFile("ClusterInfo.root","recreate");
  //hClusterE = new TH1F("clusterEnergy", "clusterEnergy", 100, 0, 100);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyCorrectionSC::MakeEnergyCorrections(const pandora::Cluster *const pCluster, float &correctedHadronicEnergy) const
{
  //std::cout << "Starts EnergyCorrectionSC" << std::endl;

    if(NULL == pCluster)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    if(0 == pCluster->GetNCaloHits())
      {
        correctedHadronicEnergy = 0.f;
        return STATUS_CODE_SUCCESS;
      }

    // PFO only has the nonIsolated calo hits associated to it.  Topologically this is right, but the energy of the
    // isolated hits is used in PFO creation.
    const pandora::OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    // This converts the OrderedCaloHitList into a standard CaloHitList
    pandora::CaloHitList nonIsolatedCaloHitList;
    orderedCaloHitList.GetCaloHitList(nonIsolatedCaloHitList);
    // Isolated calo hit list associated to the cluster
    const pandora::CaloHitList &isolatedCaloHitList(pCluster->GetIsolatedCaloHitList());
    pandora::CaloHitList clusterCaloHitList;
    clusterCaloHitList.insert(nonIsolatedCaloHitList.begin(), nonIsolatedCaloHitList.end());
    clusterCaloHitList.insert(isolatedCaloHitList.begin(), isolatedCaloHitList.end());
    
    //Cluster important information
    float clusterHADenergy = pCluster->GetHadronicEnergy();
    //const unsigned int nHitsInCluster(pCluster->GetNCaloHits());    
    //pandora::CaloHitList clusterCaloHitList;
    //pCluster->GetOrderedCaloHitList().GetCaloHitList(clusterCaloHitList);

    //if (5>= nHitsInCluster)
    //return STATUS_CODE_SUCCESS;

    //const pandora::TrackList &trackList = pCluster->GetAssociatedTrackList();
    //const unsigned int nTrackAssociations = trackList.size();

    //std::cout << "nTrackAssociations " << nTrackAssociations << std::endl;

    //if (0 == nTrackAssociations){
    //std::cout << "no track associated" << std::endl;
    //return STATUS_CODE_SUCCESS;
    //}

    //Here to check the cluster type so that the corresponding correction function is called
    bool isHCalCluster(false), isECalCluster(false);
    
    //std::cout << "The hadronic energy for this cluster is:" << pCluster->GetHadronicEnergy() << std::endl;

    //C_global distribution
    float meanHitEnergy = pCluster->GetHadronicEnergy()/clusterCaloHitList.size();
    //std::cout << "clusterEnergy " << pCluster->GetHadronicEnergy() << " ncalohits " << clusterCaloHitList.size() << " meanHitEnergy " << meanHitEnergy << std::endl;
    //std::cout << "isolated had e " << pCluster->GetIsolatedHadronicEnergy() << std::endl;

    int N_lim = 0, N_av = 0;
    float maxDiff = -1e6;
    float maxHitEnergy = -1e6;
    float minHitEnergy = 1e6;
    float sumEhit = 0;
    std::vector<float> clusterHitEnergy;
    for(pandora::CaloHitList::const_iterator iter = clusterCaloHitList.begin() , endIter = clusterCaloHitList.end() ; endIter != iter ; ++iter)
      {
	const pandora::CaloHit *pCaloHit = *iter;
	float hitEnergy = pCaloHit->GetHadronicEnergy();
	float hitEnergyinMIP = pCaloHit->GetMipEquivalentEnergy();
	clusterHitEnergy.push_back(hitEnergy);
	
	sumEhit += hitEnergy;
	
	//std::cout << "hitEnergy " << hitEnergy << " hitEnergyinMIP " << hitEnergyinMIP << std::endl;
	
	if (fabs(hitEnergy-meanHitEnergy)>maxDiff)
	  maxDiff = fabs(hitEnergy-meanHitEnergy);
	if (hitEnergy>maxHitEnergy)
	  maxHitEnergy = hitEnergy;
	if (hitEnergy<minHitEnergy)
	  minHitEnergy = hitEnergy;
	
	//std::cout << "fabs(hitEnergy-meanHitEnergy) " << fabs(hitEnergy-meanHitEnergy) << " maxDiff " << maxDiff << std::endl;
	//std::cout << "maxHitEnergy " << maxHitEnergy << std::endl;
	
	if (hitEnergy<meanHitEnergy)
	  N_av++;
	if (hitEnergyinMIP<5.)
	  N_lim++;
      }
    
    //std::cout << "N_av " << N_av << " N_lim " << N_lim << std::endl;

    float C_global(0.f);
    if (N_av==0)
      return STATUS_CODE_SUCCESS;
    else
      C_global = (float)N_lim/N_av;

    if (clusterHADenergy<3)
      return STATUS_CODE_SUCCESS;
    
    if (C_global>2 || (maxDiff/(maxHitEnergy-minHitEnergy))<0.3 )
      return STATUS_CODE_SUCCESS;

    /*
    if (clusterHADenergy<10){//Scale hot hadrons
      float m_minHitsForHotHadron(5),
	m_maxHitsForHotHadron(100),
	m_hotHadronInnerLayerCut(10),
	m_hotHadronMipFractionCut(0.4),
	m_hotHadronNHitsCut(50),
	m_hotHadronMipsPerHit(15.f),
	m_scaledHotHadronMipsPerHit(5.f);
    
      // Initial hot hadron cuts
      if ((nHitsInCluster < m_minHitsForHotHadron) || (nHitsInCluster > m_maxHitsForHotHadron))
	return STATUS_CODE_SUCCESS;
      
      if ((pCluster->GetInnerPseudoLayer() < m_hotHadronInnerLayerCut) && (pCluster->GetMipFraction() < m_hotHadronMipFractionCut) &&
	  (nHitsInCluster > m_hotHadronNHitsCut))
	{
	  return STATUS_CODE_SUCCESS;
	}
      
      // Finally, check the number of mips per hit
      float clusterMipEnergy(0.);
      const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
      
      for (OrderedCaloHitList::const_iterator layerIter = orderedCaloHitList.begin(), layerIterEnd = orderedCaloHitList.end();
	   layerIter != layerIterEnd; ++layerIter)
	{
	  for (CaloHitList::const_iterator hitIter = layerIter->second->begin(), hitIterEnd = layerIter->second->end();
	       hitIter != hitIterEnd; ++hitIter)
	    {
	      clusterMipEnergy += (*hitIter)->GetMipEquivalentEnergy();
	    }
	}
      
      const float meanMipsPerHit(clusterMipEnergy / static_cast<float>(nHitsInCluster));
      
      if ((meanMipsPerHit > 0.f) && (meanMipsPerHit > m_hotHadronMipsPerHit))
	correctedHadronicEnergy *= m_scaledHotHadronMipsPerHit / meanMipsPerHit;
      
      return STATUS_CODE_SUCCESS;
    }
    */


    this->clusterType(clusterCaloHitList, isECalCluster, isHCalCluster);

    if (isECalCluster)
    {
        //EnergyCorrectionSC::ECalClusterEnergyCorrectionFunction(clusterCaloHitList, correctedHadronicEnergy);
        return STATUS_CODE_SUCCESS;
    }

    else if (isHCalCluster){
      EnergyCorrectionSC::SCClusterEnergyCorrectionFunction(clusterHADenergy,clusterCaloHitList, correctedHadronicEnergy);
      //hClusterE->Fill(clusterHADenergy);
    }
    else
      EnergyCorrectionSC::SCHCalECalSplitClusterEnergyCorrectionFunction(clusterHADenergy,clusterCaloHitList, correctedHadronicEnergy);
    
    return STATUS_CODE_SUCCESS;

}


//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyCorrectionSC::FindDensity(const pandora::CaloHit *pCaloHit, float &hitBinDensity) const
{
  float hitEnergy = pCaloHit->GetHadronicEnergy();
  //float hitEnergyInMIP = pCaloHit->GetMipEquivalentEnergy();
  //float mip2gev = hitEnergy/hitEnergyInMIP;

  float cellvolume = (pCaloHit->GetCellSize0()) * (pCaloHit->GetCellSize1()) * (pCaloHit->GetCellThickness()) /1000000.;
  float hitEnergyDensity = hitEnergy/cellvolume;

  const int NBIN = 10;
  //float lowMIP[NBIN]  = {0.3,  2, 5.5,  8, 10, 14, 17, 21, 25, 30};
  //float highMIP[NBIN] = {  2, 5.5,   8, 10, 14, 17, 21, 25, 30, 1e6};
  float lowDensity[NBIN]  = {0, 2,   5, 7.5, 9.5, 13, 16,   20, 23.5,  28};
  float highDensity[NBIN] = {2, 5, 7.5, 9.5,  13, 16, 20, 23.5,   28,  1e6};

  for (int ibin = 0; ibin < NBIN; ibin++){
    if (hitEnergyDensity>=lowDensity[ibin] && hitEnergyDensity<highDensity[ibin]){
      hitBinDensity = (lowDensity[ibin]+highDensity[ibin])/2;
      if (ibin==(NBIN-1))
	hitBinDensity = 30;    
    }
  }

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyCorrectionSC::clusterType(const pandora::CaloHitList &caloHitList, bool &isECalCluster, bool &isHCalCluster) const
{
    int nECalHits(0), nHCalHits(0);

    for(pandora::CaloHitList::const_iterator iter = caloHitList.begin() , endIter = caloHitList.end() ; endIter != iter ; ++iter)
    {
        const pandora::CaloHit *pCaloHit = *iter;

        if (HCAL == pCaloHit->GetHitType())
            nHCalHits++;

        else if (ECAL == pCaloHit->GetHitType())
            nECalHits++;

        //std::cout << "pCaloHit->GetHitType():   " << pCaloHit->GetHitType() << std::endl;
        //std::cout << "nECalHits:                " << nECalHits << std::endl;
        //std::cout << "nHCalHits:                " << nHCalHits << std::endl;
    }

    if (nECalHits != 0 && nHCalHits == 0)
        isECalCluster = true;

    else if (nHCalHits != 0 && nECalHits == 0)
        isHCalCluster = true;

    //std::cout << "nECalHits:    " << nECalHits << std::endl;
    //std::cout << "nHCalHits:    " << nHCalHits << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyCorrectionSC::ECalClusterEnergyCorrectionFunction(const pandora::CaloHitList &caloHitList, float &energyCorrection) const
{
    std::cout << "========ECAL CONTAINED CLUSTER -> Calling EnergyCorrectionSC::ECalClusterEnergyCorrectionFunction========" << std::endl;

    for(pandora::CaloHitList::const_iterator iter = caloHitList.begin() , endIter = caloHitList.end() ; endIter != iter ; ++iter)
    {
        const pandora::CaloHit *pCaloHit = *iter;
        energyCorrection += pCaloHit->GetHadronicEnergy();
        //std::cout << "Calo Hit Energy: " << pCaloHit->GetInputEnergy() << std::endl;
    }

    //std::cout << "energyCorrection:   " << energyCorrection <<std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyCorrectionSC::SCClusterEnergyCorrectionFunction(float clusterEnergyEstimation, const pandora::CaloHitList &caloHitList, float &energyCorrection) const
{
    //std::cout << "========HCAL CONTAINED CLUSTER========" << std::endl;

    //std::cout << "m_SCEnergyConstants.at(0) : " << m_SCEnergyConstants.at(0) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(1) : " << m_SCEnergyConstants.at(1) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(2) : " << m_SCEnergyConstants.at(2) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(3) : " << m_SCEnergyConstants.at(3) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(4) : " << m_SCEnergyConstants.at(4) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(5) : " << m_SCEnergyConstants.at(5) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(6) : " << m_SCEnergyConstants.at(6) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(7) : " << m_SCEnergyConstants.at(7) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(8) : " << m_SCEnergyConstants.at(8) << std::endl;

  //std::cout << "Energy correction for neutral hadron clusters: clusterEnergyEstimation " << clusterEnergyEstimation << std::endl;
  
  float E_SC(0.f);

  float p10 = m_SCEnergyConstants1.at(0);
  float p11 = m_SCEnergyConstants1.at(1);
  float p12 = m_SCEnergyConstants1.at(2);

  float p20 = m_SCEnergyConstants1.at(3);
  float p21 = m_SCEnergyConstants1.at(4);
  float p22 = m_SCEnergyConstants1.at(5);

  float p30 = m_SCEnergyConstants1.at(6);
  float p31 = m_SCEnergyConstants1.at(7);
  float p32 = m_SCEnergyConstants1.at(8);

  //float k10 = m_SCEnergyConstants2.at(0);
  //float k11 = m_SCEnergyConstants2.at(1);
  //float k12 = m_SCEnergyConstants2.at(2);

  //float k20 = m_SCEnergyConstants2.at(3);
  //float k21 = m_SCEnergyConstants2.at(4);
  //float k22 = m_SCEnergyConstants2.at(5);

  //float k30 = m_SCEnergyConstants2.at(6);
  //float k31 = m_SCEnergyConstants2.at(7);
  //float k32 = m_SCEnergyConstants2.at(8);

  //int nHits = caloHitList.size();
  //std::cout << "Number of hits in cluster = " << nHits << std::endl;

  if (m_cheating)
    clusterEnergyEstimation = m_trueEnergy;//Cheating 

  float p1 = p10 + p11*clusterEnergyEstimation + p12*clusterEnergyEstimation*clusterEnergyEstimation;
  float p2 = p20 + p21*clusterEnergyEstimation + p22*clusterEnergyEstimation*clusterEnergyEstimation;
  float p3 = p30/(p31 + exp(p32*clusterEnergyEstimation));

  //float k1 = k10 + k11*clusterEnergyEstimation + k12*clusterEnergyEstimation*clusterEnergyEstimation;
  //float k2 = k20 + k21*clusterEnergyEstimation + k22*clusterEnergyEstimation*clusterEnergyEstimation;
  //float k3 = k30/(k31 + exp(k32*clusterEnergyEstimation));

  for(pandora::CaloHitList::const_iterator iter = caloHitList.begin() , endIter = caloHitList.end() ; endIter != iter ; ++iter)
    {
        const pandora::CaloHit *pCaloHit = *iter;
        //std::cout << "Calo Hit Energy: " << pCaloHit->GetInputEnergy() << std::endl;
        //std::cout << "Calo Hit Type:   " << pCaloHit->GetHitType() << std::endl;

	float hitEnergy = pCaloHit->GetHadronicEnergy();

	float hitBinDensity(0.f);
	this->FindDensity(pCaloHit,hitBinDensity);

	double weight;
	weight = p1*exp(p2*hitBinDensity) + p3;
	
	float hitEcorr = 0;
	hitEcorr = hitEnergy*weight;//here I have to use clusterEnergy estimation

	//std::cout << "hitEcorr " << hitEcorr << std::endl;
	E_SC += hitEcorr;
    }

    energyCorrection = E_SC;

    std::cout << "Corrected energy with SC:     " << energyCorrection <<std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyCorrectionSC::SCHCalECalSplitClusterEnergyCorrectionFunction(float clusterEnergyEstimation, const pandora::CaloHitList &caloHitList, float &energyCorrection) const
{

    //std::cout << "========SPLIT CLUSTER========" << std::endl;

    //std::cout << "m_SCEnergyConstants.at(0) : " << m_SCEnergyConstants.at(0) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(1) : " << m_SCEnergyConstants.at(1) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(2) : " << m_SCEnergyConstants.at(2) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(3) : " << m_SCEnergyConstants.at(3) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(4) : " << m_SCEnergyConstants.at(4) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(5) : " << m_SCEnergyConstants.at(5) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(6) : " << m_SCEnergyConstants.at(6) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(7) : " << m_SCEnergyConstants.at(7) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(8) : " << m_SCEnergyConstants.at(8) << std::endl;

  float E_SC(0.f);

  //std::cout << "Energy correction for split clusters: clusterEnergyEstimation " << clusterEnergyEstimation << std::endl;

  float p10 = m_SCEnergyConstants1.at(0);
  float p11 = m_SCEnergyConstants1.at(1);
  float p12 = m_SCEnergyConstants1.at(2);

  float p20 = m_SCEnergyConstants1.at(3);
  float p21 = m_SCEnergyConstants1.at(4);
  float p22 = m_SCEnergyConstants1.at(5);

  float p30 = m_SCEnergyConstants1.at(6);
  float p31 = m_SCEnergyConstants1.at(7);
  float p32 = m_SCEnergyConstants1.at(8);

  //float k10 = m_SCEnergyConstants2.at(0);
  //float k11 = m_SCEnergyConstants2.at(1);
  //float k12 = m_SCEnergyConstants2.at(2);

  //float k20 = m_SCEnergyConstants2.at(3);
  //float k21 = m_SCEnergyConstants2.at(4);
  //float k22 = m_SCEnergyConstants2.at(5);

  //float k30 = m_SCEnergyConstants2.at(6);
  //float k31 = m_SCEnergyConstants2.at(7);
  //float k32 = m_SCEnergyConstants2.at(8);

  //int nHits = caloHitList.size();
  //std::cout << "Number of hits in cluster = " << nHits << std::endl;

  if (m_cheating)
    clusterEnergyEstimation = m_trueEnergy;//Cheating 

  float p1 = p10 + p11*clusterEnergyEstimation + p12*clusterEnergyEstimation*clusterEnergyEstimation;
  float p2 = p20 + p21*clusterEnergyEstimation + p22*clusterEnergyEstimation*clusterEnergyEstimation;
  float p3 = p30/(p31 + exp(p32*clusterEnergyEstimation));

  //float k1 = k10 + k11*clusterEnergyEstimation + k12*clusterEnergyEstimation*clusterEnergyEstimation;
  //float k2 = k20 + k21*clusterEnergyEstimation + k22*clusterEnergyEstimation*clusterEnergyEstimation;
  //float k3 = k30/(k31 + exp(k32*clusterEnergyEstimation));

  for(pandora::CaloHitList::const_iterator iter = caloHitList.begin() , endIter = caloHitList.end() ; endIter != iter ; ++iter)
    {
      const pandora::CaloHit *pCaloHit = *iter;
      //std::cout << "Calo Hit Energy: " << pCaloHit->GetInputEnergy() << std::endl;
      //std::cout << "Calo Hit Type:   " << pCaloHit->GetHitType() << std::endl;

	//If HCAL: correct
        if (HCAL == pCaloHit->GetHitType())
        {
	  float hitEnergy = pCaloHit->GetHadronicEnergy();
	  float hitBinDensity(0.f);
	  this->FindDensity(pCaloHit,hitBinDensity);

	  double weight;
	  weight = p1*exp(p2*hitBinDensity) + p3;

	  //std::cout << "weigth " << weight << std::endl;
	  
	  float hitEcorr = 0;
	  hitEcorr = hitEnergy*weight;//here I have to use clusterEnergy estimation

	  E_SC += hitEcorr;
        }	
        else
        {           
	  E_SC += pCaloHit->GetHadronicEnergy();
        }
    }

    energyCorrection = E_SC;

    //std::cout << "The corrected hadronic energy for a cluster of threshold hits:" << std::endl;
    //std::cout << "E_EM = " << E_EM << " E_HAD = " << E_HAD << " E_HAD_SC " << E_HAD_SC << std::endl;
    //std::cout << "energyCorrection:     " << energyCorrection << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyCorrectionSC::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "SCEnergyConstantsHighE", m_SCEnergyConstants1));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "SCEnergyConstantsSmallE", m_SCEnergyConstants2));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "Cheating", m_cheating));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "TrueEnergy", m_trueEnergy));

    //std::cout << "====== Read Settings =====" << std::endl;
    //std::cout << "m_SCEnergyConstants.at(0) : " << m_SCEnergyConstants.at(0) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(1) : " << m_SCEnergyConstants.at(1) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(2) : " << m_SCEnergyConstants.at(2) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(3) : " << m_SCEnergyConstants.at(3) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(4) : " << m_SCEnergyConstants.at(4) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(5) : " << m_SCEnergyConstants.at(5) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(6) : " << m_SCEnergyConstants.at(6) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(7) : " << m_SCEnergyConstants.at(7) << std::endl;
    //std::cout << "m_SCEnergyConstants.at(8) : " << m_SCEnergyConstants.at(8) << std::endl;

    return STATUS_CODE_SUCCESS;
}

//I =1; J = 2, K = 12
//I = 1; j = 15; k = 12

