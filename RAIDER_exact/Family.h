#ifndef RAIDER_FAMILY
#define RAIDER_FAMILY
#include <vector>
#include "LmerVector.h"
#include <iostream>
#include <assert.h>

using namespace std;

class Family {
public:

	void adopt(LmerVector *v) {
		v->setFamily(this);
		vectors.push_back(v);
	}

	uint size() const {
		if (vectors.size() == 0) {
			return 0;
		}
		return vectors.front()->size();
	}

	uint repeatLength(uint L) const {
		return vectors.size() + L - 1;
	}

	vector<LmerVector*>* getLmers() {
		return &vectors;
	}

	LmerVector* at(uint index) {
		return vectors.at(index);
	}

	void push_back(LmerVector* v) {
		vectors.push_back(v);
	}

	LmerVector* getPrefix() const {
		return vectors.front();
	}

	LmerVector* getSuffix() const {
		return vectors.back();
	}

	vector<LmerVector*> vectors;
};

#endif //RAIDER_FAMILY
