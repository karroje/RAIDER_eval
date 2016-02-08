#ifndef RAIDER
#define RAIDER

#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <stdexcept>
#include <vector>
#include "Family.h"
#include <unordered_map>

using namespace std;

typedef unordered_map<size_t, LmerVector*> LmerMap;
typedef unordered_map<size_t, uint> LmerToIntMap;

size_t baseToInt(char b) {
	switch (toupper(b)) {
	case 'A':
		return 0;
	case 'C':
		return 1;
	case 'G':
		return 2;
	case 'T':
		return 3;
	default:
		throw std::out_of_range("Unknown base type");
		return -1;
	}
}

size_t seedToInt(const char* seed, uint L) {
	//TODO use bitset or vector<bool> to allow for potentially large bitsets
	size_t result = 0;
       	for (uint i = 0; i < L; i++) {
	  result <<= 2;
	  result |= baseToInt(seed[i]);
	}
	return result;
}

bool getNextSeed(seqan::Dna5String &sequence, uint &index, uint seqLength, uint L, char* seed) {
	seed[L - 1] = 0;
	for (uint i = 0; i < L && index + i < seqLength; i++) {
		if (toupper((char) sequence[index + i]) == 'N') {
			index = index + i + 1;
			i = -1;
		} else {
			seed[i] = sequence[index + i];
		}
	}
	return seed[L - 1];
}

bool isNewFamily(Family* fam, LmerVector* v) {
	assert(!v->getFamily());
	if (!fam || fam->size() > 2) {
		return true;
	}
	assert(fam->size() == 2);
	LmerVector* suffix = fam->getSuffix();
	uint s1 = suffix->front();
	uint s2 = suffix->back();
	uint prev = v->front();
	uint index = v->back();
	// Ensure both Lmer pairs are the same distance apart
	if (index - prev != s2 - s1) {
		return true;
	}
	// Ensure first Lmer pair overlaps second Lmer pair
	if (index != s2 + 1) {
		return true;
	}
	return false;
}


Family* splitRepeatsByLmer(Family* fam, LmerVector *v, bool keepV) {
	assert(fam == v->getFamily());
	vector<LmerVector*> *lmers = fam->getLmers();
	Family *newFam = new Family;
	uint oldLength = lmers->size();
	uint offset = oldLength;
	for (uint i = 0; i < oldLength; i++) {
		LmerVector *u = lmers->at(i);
		if (u == v) {
			offset = i;
		}
		if ((!keepV && i >= offset) || i > offset) {
			newFam->adopt(u);
		}
	}
	lmers->resize(keepV ? offset + 1 : offset);
	lmers->shrink_to_fit();
	return newFam;
}

bool isPrefix(LmerVector *v) {
	return v->getFamily()->getPrefix() == v;
}

LmerVector* getAndInsert(LmerMap &lmers, const size_t seed, const uint index) {
	auto it = lmers.find(seed);
	if (it == lmers.end()) {
		return nullptr;
	}
	LmerVector *v = (*it).second;
	v->push_back(index);
	return v;
}

void setLmers(seqan::Dna5String &sequence, const uint seqLength, const uint MIN, const uint L, LmerMap &lmers) {
	LmerToIntMap counts;
	uint index = -1;
	char seed[L];

	while (getNextSeed(sequence, ++index, seqLength, L, seed)) {
		size_t key = seedToInt(seed, L);
		counts[key]++;
	}
	for (auto kv : counts) {
	  if (kv.second >= MIN) {
	    lmers[kv.first] = new LmerVector();
	  }
	}
}

void getElementaryFamilies(seqan::Dna5String &sequence, const uint MIN, const uint L, vector<Family*> &families) {
	const uint seqLength = seqan::length(sequence);

	uint index = -1;
	char seed[L];
	LmerMap lmers;
	Family* fam = nullptr;
	Family* newFam = nullptr;
	uint c = 0;

	setLmers(sequence, seqLength, MIN, L, lmers);

	while (getNextSeed(sequence, ++index, seqLength, L, seed)) {
		size_t key = seedToInt(seed, L);
		LmerVector *v = getAndInsert(lmers, key, index);
		if (v == nullptr) {
			continue;
		}
		if (fam) {
			if (v != fam->at(c)) {
				families.push_back(splitRepeatsByLmer(fam, fam->at(c), false));
				fam = nullptr;
			} else {
				c++;
			}
		}
		if (v->size() == 2) {
			if (isNewFamily(newFam, v)) {
				newFam = new Family();
				families.push_back(newFam);
			}
			newFam->adopt(v);
		} else if (v->size() > 2) {
			Family* vFam = v->getFamily();
			if (v == vFam->getPrefix()) {
				fam = vFam;
				c = 1;
			} else if (fam != vFam) {
				fam = splitRepeatsByLmer(vFam, v, false);
				families.push_back(fam);
				c = 1;
			}
			if (fam && fam->getSuffix() == v) {
				fam = nullptr;
			}
		}
	}
}

#endif //RAIDER
