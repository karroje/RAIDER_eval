#ifndef SCANER_FAMILY
#define SCANER_FAMILY
#include <vector>
#include "LmerVector.h"
#include <iostream>
//#include <assert.h>

using namespace std;

class Family {
public:
	Family(){

	}

	Family(uint L, LmerVector *v){
		repLength = L;
		adopt(v);
	}

	// A family adopts an LmerVector. The LmerVector is set with family
	// The family adds LmerVector to end of list of vectors.
	// Expected end is updated as last location of v + 1
	void adopt(LmerVector *v) {
		v->setFamily(this);
		vectors.push_back(v);
		setLast(v);
		setExpectedEnd(v->back() + 1);
	}

	uint size() const {
		if (vectors.size() == 0) {
			return 0;
		}
		return vectors.front()->size();
	}

	// no longer true necessarily (can skip lmers)
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

	LmerVector* getLast() {
		return last;
	}

	uint getLastIndex() {
		return last_index;
	}

	void setLast(LmerVector* v) {
		assert(v->getFamily() == this);
		last = v;
		last_index = v->back();
	}

	bool lastRepeatComplete() const {
		return last == getSuffix();
	}

	void setExpectedEnd(uint expected) {
		expected_end = expected;
	}

	uint getExpectedEnd() {
		return expected_end;
	}

	uint extendRepeatLength(uint offset){
		repLength += offset;
		return repLength;
	}

	uint getRepeatLength(){
		return repLength;
	}

	bool hasConflict(LmerVector* v, uint L){
		uint i = 0;
		uint j = 0;
		while(i < last->size() && j < v->size()){
			// start of v location after end of last
			if((*v)[j] > (*last)[i] + L){
				i++;
			}
			// start of last after end of v
			else if((*last)[i] - repLength > (*v)[j] + L){
				j++;
			}
			else{
				return true;
			}
		}
		return false;
	}

	uint repLength;
	LmerVector* last;
	uint last_index;
	uint expected_end;
	vector<LmerVector*> vectors;
};

#endif //SCANER_FAMILY
