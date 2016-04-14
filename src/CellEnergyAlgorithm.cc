/**
 *  @file   MarlinPandora/src/CellEnergyAlgorithm.cc
 * 
 *  @brief  Implementation of the cell energy algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "PandoraMonitoringApi.h"

#include "CellEnergyAlgorithm.h"

using namespace pandora;

//------------------------------------------------------------------------------------------------------------------------------------------

CellEnergyAlgorithm::CellEnergyAlgorithm() :
    m_RootFileName("Test.root")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

CellEnergyAlgorithm::~CellEnergyAlgorithm()
{   
    std::cout << "Exiting Cell Energy Algorithm and saving root file : " << m_RootFileName << std::endl;
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ClusterEnergyTree", m_RootFileName, "UPDATE"));
    //PandoraMonitoringApi::SaveTree(this->GetPandora(), )
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CellEnergyAlgorithm::Run()
{
    std::cout << "===== Running Cell Energy Algorithm =====" << std::endl;

    // Algorithm code here
    const pandora::PfoList *pPfoList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pPfoList));

    std::vector<float> *pfoMomentum = new std::vector<float>;
    std::vector<float> *pfoCosTheta = new std::vector<float>;

    for (pandora::PfoList::const_iterator pfoIter = pPfoList->begin(), pfoIterEnd = pPfoList->end(); pfoIter != pfoIterEnd; ++pfoIter)
    {
        const pandora::Pfo *const pPfo = *pfoIter;
        const CartesianVector cartesianVector = pPfo->GetMomentum();
        const float px = cartesianVector.GetX();
        const float py = cartesianVector.GetY();
        const float pz = cartesianVector.GetZ();
        //std::cout << px << " " << py << " " << pz << std::endl; 
        const float momentum(std::sqrt(px * px + py * py + pz * pz));
        const float cosTheta((momentum > std::numeric_limits<float>::epsilon()) ? pz / momentum : -999.f);
        pfoMomentum->push_back(momentum);
        pfoCosTheta->push_back(cosTheta);
    }

    int numberOfPfos(pPfoList->size());

    const pandora::ClusterList *pCurrentClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCurrentClusterList));

    std::vector<int> *numberOfHitsInCluster = new std::vector<int>;
    std::vector<float> *rawEnergyOfCluster = new std::vector<float>;
    std::vector<float> *rawIsolatedEnergyOfCluster = new std::vector<float>;
    std::vector<float> *correctedEnergyOfCluster = new std::vector<float>;
    std::vector<int> *nECalHits = new std::vector<int>;
    std::vector<int> *nHCalHits = new std::vector<int>;
    std::vector<int> *isEMShower = new std::vector<int>;
    std::vector<float> *Cglobal = new std::vector<float>;
    std::vector<float> *MeanHitEnergy = new std::vector<float>;
    std::vector<float> *MaxDiff = new std::vector<float>;
    std::vector<float> *MaxHitEnergy = new std::vector<float>;
    std::vector<float> *MinHitEnergy = new std::vector<float>;
    //std::vector<std::vector<float> > *ClusterHitEnergy = new std::vector<std::vector<float> >;

    for (pandora::ClusterList::const_iterator clusterIter = pCurrentClusterList->begin(), clusterIterEnd = pCurrentClusterList->end(); clusterIter != clusterIterEnd; ++clusterIter)
    {
        const Cluster *const pCluster = *clusterIter;
        const bool emShower(PandoraContentApi::GetPlugins(*this)->GetParticleId()->IsEmShower(pCluster));
        int necalHits(0);
        int nhcalHits(0);
        //const float clusterHadEnergy = pCluster->GetHadronicEnergy();

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
    
        this->ClusterType(clusterCaloHitList,necalHits,nhcalHits);

//        if (!isEMShower)
        std::cout << "Number of hits in cluster : " << pCluster->GetNCaloHits() << std::endl;
        numberOfHitsInCluster->push_back(clusterCaloHitList.size());
        rawEnergyOfCluster->push_back(pCluster->GetHadronicEnergy());
        correctedEnergyOfCluster->push_back(pCluster->GetCorrectedHadronicEnergy(this->GetPandora()));
        nECalHits->push_back(necalHits);
        nHCalHits->push_back(nhcalHits);
        isEMShower->push_back(emshower);
        rawIsolatedEnergyOfCluster->push_back(pCluster->GetIsolatedHadronicEnergy());

	//C_global distribution
	float meanHitEnergy = pCluster->GetHadronicEnergy()/clusterCaloHitList.size();
	//std::cout << "clusterEnergy " << pCluster->GetHadronicEnergy() << " ncalohits " << clusterCaloHitList.size() << " meanHitEnergy " << meanHitEnergy << std::endl;
	//std::cout << "isolated had e " << pCluster->GetIsolatedHadronicEnergy() << std::endl;

	int N_lim = 0, N_av = 0;
	float maxDiff = -1e6;
	float maxHitEnergy = -1e6;
	float minHitEnergy = 1e6;
	float sumEhit = 0;
	//std::vector<float> *clusterHitEnergy = new std::vector<float>;
	for(pandora::CaloHitList::const_iterator iter = clusterCaloHitList.begin() , endIter = clusterCaloHitList.end() ; endIter != iter ; ++iter)
	  {
	    const pandora::CaloHit *pCaloHit = *iter;
	    float hitEnergy = pCaloHit->GetHadronicEnergy();
	    float hitEnergyinMIP = pCaloHit->GetMipEquivalentEnergy();
	    //clusterHitEnergy->push_back(hitEnergy);

	    sumEhit += hitEnergy;

	    std::cout << "hitEnergy " << hitEnergy << " hitEnergyinMIP " << hitEnergyinMIP << std::endl;
	    
	    if (fabs(hitEnergy-meanHitEnergy)>maxDiff)
	      maxDiff = fabs(hitEnergy-meanHitEnergy);
	    if (hitEnergy>maxHitEnergy)
	      maxHitEnergy = hitEnergy;
	    if (hitEnergy<minHitEnergy)
	      minHitEnergy = hitEnergy;

	    std::cout << "fabs(hitEnergy-meanHitEnergy) " << fabs(hitEnergy-meanHitEnergy) << " maxDiff " << maxDiff << std::endl;
	    std::cout << "maxHitEnergy " << maxHitEnergy << std::endl;

	    if (hitEnergy<meanHitEnergy)
	      N_av++;
	    if (hitEnergyinMIP<5.)
	      N_lim++;
	  }
	
	std::cout << "sumEhit " << sumEhit << std::endl;
	
	if ((float)N_lim/N_av==1)
	  std::cout << "N_av " << N_av << " N_lim " << N_lim << " N_lim/N_av " << (float)N_lim/N_av << std::endl;

	if (maxDiff/maxHitEnergy>0.8){
	  std::cout << "maxDiff/maxHitEnergy>0.8 " << std::endl;
	  std::cout << "maxDiff: " << maxDiff << " maxHitEnergy: " << maxHitEnergy << " minHitEnergy " << minHitEnergy << std::endl;
	  std::cout << "maxDiff in MIP: " << maxDiff/0.0225 << " maxHitEnergy in MIP: " << maxHitEnergy/0.0225 << std::endl;
	} else {
	  std::cout << "maxDiff/maxHitEnergy<0.8 " << std::endl;
	  std::cout << "maxDiff: " << maxDiff << " maxHitEnergy: " << maxHitEnergy << " minHitEnergy " << minHitEnergy << std::endl;
	  std::cout << "maxDiff in MIP: " << maxDiff/0.0225 << " maxHitEnergy in MIP: " << maxHitEnergy/0.0225 << std::endl;
	}

	if (N_av>0)
	  Cglobal->push_back((float)N_lim/N_av);	
	else
	  Cglobal->push_back(1e6);
	MeanHitEnergy->push_back(meanHitEnergy);
	MaxDiff->push_back(maxDiff);
	MaxHitEnergy->push_back(maxHitEnergy);
	MinHitEnergy->push_back(minHitEnergy);
	//ClusterHitEnergy->push_back(clusterHitEnergy);

        std::string detectorView = "default";
        //PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, (detectorView.find("xz") != std::string::npos) ? DETECTOR_VIEW_XZ : (detectorView.find("xy") != std::string::npos) ? DETECTOR_VIEW_XY : DETECTOR_VIEW_DEFAULT, -1.f, clusterHadEnergy);
        //PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &clusterCaloHitList, "SoftCompCaloHits", SOFTCOMPWEIGHT);
        //PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &clusterCaloHitList, "SoftCompCaloHits", CLUSTERHADE);
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "RawEnergyOfCluster", rawEnergyOfCluster));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "RawIsolatedEnergyOfCluster", rawIsolatedEnergyOfCluster));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "CorrectedEnergyOfCluster", correctedEnergyOfCluster));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "NumberOfHitsInCluster", numberOfHitsInCluster));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "NumberOfPFOsInEvent", numberOfPfos));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "PFOMomentum", pfoMomentum));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "PFOCosTheta", pfoCosTheta));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "nECalHits", nECalHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "nHCalHits", nHCalHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "IsEMShower", isEMShower));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "Cglobal", Cglobal));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "MeanHitEnergy", MeanHitEnergy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "MaxDiff", MaxDiff));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "MaxHitEnergy", MaxHitEnergy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "MinHitEnergy", MinHitEnergy)); 
    //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ClusterEnergyTree", "ClusterHitEnergy", &ClusterHitEnergy));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "ClusterEnergyTree"));

   return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CellEnergyAlgorithm::ClusterType(const pandora::CaloHitList &caloHitList, int &nECalHits, int &nHCalHits) const
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

/*StatusCode CellEnergyAlgorithm::Display(const CaloHitList &fullCaloHitList, const CaloHitList &hotCaloHitList) const
{
    std::string detectorView = "default";
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, (detectorView.find("xz") != std::string::npos) ? DETECTOR_VIEW_XZ : (detectorView.find("xy") != std::string::npos) ? DETECTOR_VIEW_XY : DETECTOR_VIEW_DEFAULT, -1.f, 1.f);

//    PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), pPfoList, "AllPFOs", BLACK, true, true);

    PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &fullCaloHitList, "AllCaloHits", GREEN);
    PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &hotCaloHitList, "HotCaloHits", RED);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    return STATUS_CODE_SUCCESS;
}*/

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CellEnergyAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    // Read settings from xml file here

//    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "HotHadCellEnergy", m_HotHadCellEnergy));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_RootFileName));

    return STATUS_CODE_SUCCESS;
}
