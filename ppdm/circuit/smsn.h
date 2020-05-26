#ifndef __SMSN_CIRCUIT_H__BY_YSSHIN_
#define __SMSN_CIRCUIT_H__BY_YSSHIN_

#include "circuit.h"

class CSmsn : public CCircuit
{
public:
	BOOL Create(int nParties, const vector<int>& params);

public:
	void PutAdditionLayer();
	void PutOutputLayer();

private:
	int m_nInputs;
	int m_nRep;
	vector<int>	m_vLayerInputs;
	vector<int> m_vLayerCompares;
	vector<int> m_vLayerOutputs;
};

#endif // __SMSN_CIRCUIT_H__BY_YSSHIN_

