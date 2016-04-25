/**
 *  @file   MarlinPandora/src/RegisterHitsForSC.cc
 * 
 *  @brief  Implementation of the register hits for SC algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "PandoraMonitoringApi.h"

#include "RegisterHitsForSC.h"

using namespace pandora;

//------------------------------------------------------------------------------------------------------------------------------------------

RegisterHitsForSC::RegisterHitsForSC() :
    m_myRootFileName("HitsForSC.root")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

RegisterHitsForSC::~RegisterHitsForSC()
{   
    std::cout << "Exiting RegisterHitsForSC algorithm and saving root file : " << m_myRootFileName << std::endl;
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "HitEnergyTree", m_myRootFileName, "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RegisterHitsForSC::Run()
{
    std::cout << "===== Running RegisterHitsForSC algorithm =====" << std::endl;

    // Add PFO information
    //********************
    const pandora::PfoList *pPfoList = NULL;
    std::vector<float> *EnergyOfPfos = new std::vector<float>;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pPfoList));
    for (pandora::PfoList::const_iterator pfoIter = pPfoList->begin(), pfoIterEnd = pPfoList->end(); pfoIter != pfoIterEnd; ++pfoIter)
      {
	const pandora::Pfo *const pPfo = *pfoIter;
	EnergyOfPfos->push_back(pPfo->GetEnergy());
      }

    int numberOfPfos(pPfoList->size());
    //********************

    std::cout << "numberOfPfos " << numberOfPfos << std::endl;

    //Add cluster and information
    //********************
    const pandora::ClusterList *pCurrentClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCurrentClusterList));

    int numberOfClusters(pCurrentClusterList->size());
    std::vector<int> *numberOfHitsInCluster = new std::vector<int>;
    std::vector<float> *rawEnergyOfCluster = new std::vector<float>;
    std::vector<float> *rawIsolatedEnergyOfCluster = new std::vector<float>;
    std::vector<int> *isEMShower = new std::vector<int>;
    std::vector<int> *nECalHits = new std::vector<int>;
    std::vector<int> *nHCalHits = new std::vector<int>;
    std::vector<float> *m_CellSize0 = new std::vector<float>;
    std::vector<float> *m_CellSize1 = new std::vector<float>;
    std::vector<float> *m_CellThickness = new std::vector<float>;
    std::vector<float> *HitEnergies = new std::vector<float>;
    std::vector<int> *hitType = new std::vector<int>;

    for (pandora::ClusterList::const_iterator clusterIter = pCurrentClusterList->begin(), clusterIterEnd = pCurrentClusterList->end(); clusterIter != clusterIterEnd; ++clusterIter)
    {
        const Cluster *const pCluster = *clusterIter;
        const bool emShower(PandoraContentApi::GetPlugins(*this)->GetParticleId()->IsEmShower(pCluster));
        int necalHits(0);
        int nhcalHits(0);

        int emshower(0);
        if (emShower) emshower = 1;

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
    
	for(pandora::CaloHitList::const_iterator hitIter = clusterCaloHitList.begin() , endhitIter = clusterCaloHitList.end() ; endhitIter != hitIter ; ++hitIter)
	  {   
	    const pandora::CaloHit *pCaloHit = *hitIter;
	    m_CellSize0->push_back(pCaloHit->GetCellSize0());
	    m_CellSize1->push_back(pCaloHit->GetCellSize1());
	    m_CellThickness->push_back(pCaloHit->GetCellThickness());


	    HitEnergies->push_back(pCaloHit->GetHadronicEnergy());
	    if (HCAL == pCaloHit->GetHitType())
	      {
		hitType->push_back(2);
	      }
	    
	    else if (ECAL == pCaloHit->GetHitType())
	      {
		hitType->push_back(1);
	      }
	    
	    else
	      {
		hitType->push_back(3);
	      }
	  }
	
	
        this->ClusterType(clusterCaloHitList,necalHits,nhcalHits);

        std::cout << "Number of hits in cluster : " << clusterCaloHitList.size() << std::endl;
	std::cout << "number of nonIsolated calo hits " << nonIsolatedCaloHitList.size() << std::endl;
	std::cout << "number of isolated calo hits " << isolatedCaloHitList.size() << std::endl;
	std::cout << "necalHits " << necalHits << " nhcalHits " << nhcalHits << std::endl;

	if (clusterCaloHitList.size() != (necalHits+nhcalHits))
	  std::cout << "clusterCaloHitList.size() != (necalHits+nhcalHits) " << (necalHits+nhcalHits) << std::endl;

	//std::vector<float> ecalHitEnergies, hcalHitEnergies;
	//this->ExtractCaloHits(clusterCaloHitList, ecalHitEnergies, hcalHitEnergies);	

        numberOfHitsInCluster->push_back(clusterCaloHitList.size());
        rawEnergyOfCluster->push_back(pCluster->GetHadronicEnergy());
        rawIsolatedEnergyOfCluster->push_back(pCluster->GetIsolatedHadronicEnergy());
        isEMShower->push_back(emshower);
        nECalHits->push_back(necalHits);
        nHCalHits->push_back(nhcalHits);
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "numberOfPfos", numberOfPfos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "numberOfClusters", numberOfClusters));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "NumberOfHitsInCluster", numberOfHitsInCluster));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "EnergyOfPfos", EnergyOfPfos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "RawEnergyOfCluster", rawEnergyOfCluster));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "RawIsolatedEnergyOfCluster", rawIsolatedEnergyOfCluster));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "IsEMShower", isEMShower)); 
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "nECalHits", nECalHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "nHCalHits", nHCalHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "CellSize0", m_CellSize0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "CellSize1", m_CellSize1));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "CellThickness", m_CellThickness));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "HitEnergies", HitEnergies));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "HitType", hitType));
    //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "ECalHitEnergies", ECalHitEnergies));
    //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "HitEnergyTree", "HCalHitEnergies", HCalHitEnergies));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "HitEnergyTree"));

   return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RegisterHitsForSC::ClusterType(const pandora::CaloHitList &caloHitList, int &nECalHits, int &nHCalHits) const
{
    for(pandora::CaloHitList::const_iterator iter = caloHitList.begin() , endIter = caloHitList.end() ; endIter != iter ; ++iter)
    {   
        const pandora::CaloHit *pCaloHit = *iter;
        
        if (HCAL == pCaloHit->GetHitType())
            nHCalHits++;
        
        else if (ECAL == pCaloHit->GetHitType())
            nECalHits++;
    }
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RegisterHitsForSC::ExtractCaloHits(const pandora::CaloHitList &caloHitList, std::vector<float> &ECalHitEnergies, std::vector<float> &HCalHitEnergies) const
{
    for(pandora::CaloHitList::const_iterator iter = caloHitList.begin() , endIter = caloHitList.end() ; endIter != iter ; ++iter)
    {   
        const pandora::CaloHit *pCaloHit = *iter;
	const float hitEnergy(pCaloHit->GetHadronicEnergy());
        
        if (HCAL == pCaloHit->GetHitType())
	  HCalHitEnergies.push_back(hitEnergy);
        
        else if (ECAL == pCaloHit->GetHitType())
	  ECalHitEnergies.push_back(hitEnergy);
    }
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RegisterHitsForSC::ReadSettings(const TiXmlHandle xmlHandle)
{
    // Read settings from xml file here

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MyRootFileName", m_myRootFileName));

    return STATUS_CODE_SUCCESS;
}
