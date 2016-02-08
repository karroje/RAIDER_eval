#ifndef SCANER_FAMILY
#define SCANER_FAMILY
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
		setLast(v);
		setExpectedEnd(v->back() + 1);
	}

    void adopt(LmerVector *v, uint L){
        v->setFamily(this);
        v->setPosition(vectors.size());
        //cout << L << endl;
        if (vectors.size() > 0){
            uint off = v->front() - vectors.back()->front();
            addToOffset(off);
            setRepeatLength(getRepeatLength() + off);
        }
        else{
            setOffset(0);
            setRepeatLength(L);
            v->setPosition(0);
        }
        vectors.push_back(v);
        setLast(v);
        setExpectedEnd(v->back() + L);
    
    }


	uint size() const {
		if (vectors.size() == 0) {
			return 0;
		}
		return vectors.front()->size();
	}

    void addToOffset(uint off){
        offset += off;
    }

    void setOffset(uint off){
        offset = off;
    }

    uint getOffset() const {
        return offset;
    }

    uint getRepeatLength() const {
        return repeatlength;
    }

    void setRepeatLength(uint L){
        repeatlength = L;
    }

	uint repeatLength(uint L) const {
		return repeatlength + L - L;
        //return vectors.size() + L - 1;
	}

	vector<LmerVector*>* getLmers() {
		return &vectors;
	}

	LmerVector* at(uint index) {
		return vectors.at(index);
	}

    LmerVector* see(uint index) const {
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

	uint getLastIndex() const {
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

	uint getExpectedEnd() const {
		return expected_end;
	}

    friend ostream &operator<<( ostream &output, const Family &F) { 
        output << "Size: " << F.vectors.size() << " Repeat Length: " << F.repeatlength << " Expected End: " << F.expected_end << " Offset: " << F.offset <<  endl;
        output << "\t\tLast:\t" << *(F.last) << endl;
        output << "\t\tPrefix:\t" << *(F.getPrefix()) << endl;
        output << "\t\tSuffix:\t" << *(F.getSuffix()) << endl;
        for(uint i = 0; i < F.vectors.size(); i++){ //F.LmerVector* l : F.getLmers()){
            output << "\t\t\t"<< i << ":\t" << *(F.see(i)) << endl;
        }
        return output;            
    }

	LmerVector* last;
	uint last_index;
	uint expected_end;
	uint repeatlength;
    uint offset;
    vector<LmerVector*> vectors;
};



#endif //SCANER_FAMILY
