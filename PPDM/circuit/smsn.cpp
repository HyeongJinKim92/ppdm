// snmn.cpp

#include "smsn.h"
#include <cassert>
using namespace std;

BOOL CSmsn::Create(int nParties, const vector<int>& vParams)
{
	// m_nNumParties = nParties;

	// parse params
	if (vParams.size() < (unsigned)2)
	{
		cout << "Error! This circuit needs " << 2
			<< "parameters: InputSize BitsNum"
			<< endl;
		return FALSE;
	}

	m_nInputs = vParams[0];
	m_nRep = vParams[1];
	m_vNumVBits.resize(m_nInputs, m_nRep); // Input수만큼 BitsNum을 할당함

	// gates for inputs
	m_vInputStart.resize(m_nInputs);
	m_vInputEnd.resize(m_nInputs);

	m_nFrontier = 2;
	for (int i = 0; i < m_nInputs; i++) // Input Gate 셋팅
	{
		m_vInputStart[i] = m_nFrontier;
		m_nFrontier += m_nRep;
		m_vInputEnd[i] = m_nFrontier - 1;
	}

	//===========================================================================
	// computes roughly the number of gates beforehand --- need this for allocating space
	// gates for each equality gate:
	// bitwise xor, bitwise negation, AND of each bitwise outputs
	int gates_st = 4 * m_nRep;
	int gates_add = 6 * m_nRep; // 6 * 16 = 96
	int gates_noninput = m_nInputs * gates_add + (m_nInputs / 2) * gates_st + m_nInputs * m_nRep; // 4 * 96 + 4 * 16

	m_othStartGate = m_nFrontier;
	m_nNumGates = m_nFrontier + gates_noninput;
	m_pGates = new GATE[m_nNumGates];
	m_nNumXORs = 0;

	//cout << "oth Start Gate " << m_othStartGate << endl;
	GATE* gate;

	for (int i = 0; i<2; i++) // GATE 1, 2 Reserved
	{
		gate = m_pGates + i;
		gate->p_ids = 0;
		gate->p_num = 0;
		gate->left = -1;
		gate->right = -1;
		gate->type = 0;
	}

	for (int i = 2; i<m_othStartGate; i++)
	{
		gate = m_pGates + i;
		gate->p_ids = 0;
		gate->p_num = 0;
		gate->left = -1;
		gate->right = -1;
		gate->type = 16;
	}

	// Prepare 0-th input layer: 
	m_vLayerInputs.resize(m_nInputs, 2); // 2 GATE IS INPUT GATE
	m_vLayerCompares.resize(m_nInputs / 2, 2); // Compare는 입력의 절반 크기
	m_vLayerOutputs.resize(1, 2); // Output은 한번만

	m_vLayerInputs = m_vInputStart;

	for (int i = 0; i < m_nInputs / 2; i++) // INPUT의 절반 크기 만큼 ADD 게이트 설정
	{
		// cout << i << "th offsets (" << m_vLayerInputs[i * 2] << " , " << m_vLayerInputs[i * 2 + 1] << ")" << endl;
		// PutAdditionLayer();
		m_vLayerCompares[i] = PutAddGate(m_vLayerInputs[i * 2], m_vLayerInputs[i * 2 + 1], m_nRep);
	}

	m_vLayerOutputs[0] = PutGTGate(m_vLayerCompares[0], m_vLayerCompares[1], m_nRep);

	PutOutputLayer();
	PatchParentFields();
	return TRUE;
};

void CSmsn::PutOutputLayer()
{

	for (int i = 0; i<1; i++)
	{
		int o_start = m_nFrontier;
		for (int j = 0; j<m_nRep; j++)
		{
			PutXORGate(m_vLayerOutputs[i] + j, 0);
		}
		int o_end = m_nFrontier - 1;
		m_vOutputStart.resize(i + 1, o_start);
		m_vOutputEnd.resize(i + 1, o_end);
	}
	m_nNumGates = m_nFrontier;
}