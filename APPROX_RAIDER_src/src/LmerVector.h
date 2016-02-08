#ifndef LMER_VECTOR
#define LMER_VECTOR
#include "LmerVector.h"

using namespace std;

#ifndef uint
typedef unsigned int uint;
#endif

class Family;

class LmerVector {
public:

	LmerVector(uint mask) {
		family = nullptr;
		seed = mask;
	}

	LmerVector(uint mask, uint position) {
		family = nullptr;
		seed = mask;
		push_back(position);
	}

	LmerVector(uint mask, vector<uint> positions) {
		family = nullptr;
		seed = mask;
        //cout << "New Lmer Vector\t";
		for(vector<uint>::iterator it = positions.begin(); it != positions.end(); ++it) {
			push_back(*it);
            //cout << (*it) << "\t";
		}
        //cout << endl << endl;
	}

	void push_back(uint val) {
		lmers.push_back(val);
	}

	uint operator [](uint i) const {
		return lmers[i];
	}

	uint & operator [](uint i) {
		return lmers[i];
	}

	uint back() {
		return lmers.back();
	}

	uint front() const {
		return lmers.front();
	}

	uint size() const {
		return lmers.size();
	}

	Family* getFamily() const {
		return family;
	}

	void setFamily(Family *fam) {
		family = fam;
	}

	uint seed;
	vector<uint> lmers;
	Family* family;
};

#endif //LMER_VECTOR
