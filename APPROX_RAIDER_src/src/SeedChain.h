#ifndef SCANER_SEEDCHAIN
#define SCANER_SEEDCHAIN

#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <stdexcept>
#include <vector>
#include "Family.h"
#include "CompareLmerVector.h"
#include <unordered_map>
#include <queue>
#include <set>

using namespace std;

typedef unordered_map<size_t, LmerVector*> LmerMap;
typedef unordered_map<size_t, vector<uint>> LmerToLocationsMap;
typedef set<size_t> SeedSet;
typedef priority_queue<LmerVector*, vector<LmerVector*>, CompareLmerVector> LmerHeap;

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

size_t seedToInt(const char* seed, const seqan::CharString &mask) {
    //TODO use bitset or vector<bool> to allower for potentially large bitsets
    size_t result = 0;
    uint L = seqan::length(mask);
    for (uint i = 0; i < L; i++) {
        if (mask[i] == '1') {
            result <<= 2;
            result |= baseToInt(seed[i]);
        }
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


uint maxLength(const vector<seqan::CharString> &masks) {
    uint max = 0;
    for (uint i = 0; i < masks.size(); i++) {
        uint len = seqan::length(masks[i]);
        if (len > max) {
            max = len;
        }
    }
    return max;
}


void setLmers(seqan::Dna5String &sequence, const uint seqLength, const seqan::CharString &mask, const uint MIN, const uint L, LmerHeap &lmers) {
    LmerToLocationsMap counts; // map seed to start and count
    SeedSet seeds; // set containing all seeds in lmers so far
    char unmasked[L]; // empty char[L]
    uint p = 3;

    for (uint i = 0; i < p; i++) {
        uint sectorMax = (i+1)*(seqLength/p);
        for (uint index = i*(seqLength/p); index < seqLength; index++) {
            if (getNextSeed(sequence, index, seqLength, L, unmasked)) {
                // unmasked now contains lmer starting at index
                size_t seed = seedToInt(unmasked, mask);
                // CS: added so don't consider previously counted
                if (seeds.count(seed) > 0){
                    continue;
                }
                if (index < sectorMax || counts.count(seed) > 0) {
                    //if(counts.count(seed) == 0){
                        //counts[seed] = new vector<int>();//make_pair(index, 1);
                    //}
                    counts[seed].push_back(index);
                }
            }
        }

        for (auto kv : counts) {
            // if frequency count is greater than min and not already in lmers
            if (kv.second.size() >= MIN){
                //LmerVector *v = new LmerVector(kv.first,kv.second);
                //lmers.push(v);
                lmers.push(new LmerVector(kv.first, kv.second)); //[kv.first] = new LmerVector();
            //    seeds.insert(kv.first);
            }
        }
        counts.clear();
    }
    
}

bool hasConflict(vector<Family*> &families, const uint L, LmerVector *u){
    for(vector<Family*>::iterator it = families.begin(); it != families.end(); ++it) {
        if((*it)-> hasConflict(u,L)){
            return true;
        }
    }
    return false;
}

bool extend(Family* fam, LmerVector* u, const uint L){
    LmerVector* v = fam->getLast();
    
    uint offset = u->front() - v->front();
    //cout <<"Offset: "<< offset << endl;
    if(offset > L){
        return false;
    }
    for(uint i = 1; i < v->size(); i++){
        if((*u)[i] - (*v)[i] != offset){
            return false;
        }
    }
    fam->adopt(u);
    fam->extendRepeatLength(offset);
    cout <<"Extending Lmer Length by "<< offset << endl;
    return true;
}



void createNewFamily(LmerHeap &lmers, const uint L, vector<Family*> &families){
    LmerVector* v = lmers.top();
    lmers.pop();
    // if(hasConflict(families,L,v)){
    //     return;
    // }
    uint currFreq = v->size();
    //uint currLoc = v->front();
    Family* fam = new Family(L,v);
    vector<LmerVector*> toReinsert;
    // while next lmer occurs with same frequency and is within L of current lmer
    while(!lmers.empty() && lmers.top()->size() == currFreq && lmers.top()->front() < fam->getLast()->front() + L){
        // cout << lmers.top()->front() << "\t" << fam->getLast()->front() << endl;
        LmerVector* u = lmers.top();
        lmers.pop();
        // if(hasConflict(families,L,u)){
        //     break;
        // }
        // Check to see if u can be used to extend length of current family
        if(!extend(fam, u, L)){
            toReinsert.push_back(u);
            // currLoc = fam->getLast()->front();
        }
        //else{
        //    toReinsert.push_back(u);
        //}
        //lmers.pop();
        //currLoc = fam->getLast()->front();
    }
    for(vector<LmerVector*>::iterator it = toReinsert.begin(); it != toReinsert.end(); ++it) {
        lmers.push(*it);
    }
    families.push_back(fam);
}





void getElementaryFamilies(seqan::Dna5String &sequence, vector<seqan::CharString> &masks, const uint MIN, vector<Family*> &families) {
    const uint seqLength = seqan::length(sequence);
    // this is for use of masks where only 1 element
    const seqan::CharString mask = masks.front();
    // default min length for lmer is length of mask
    const uint L = seqan::length(mask);
    // finds max length of masks (now irrelevant bc only 1 mask)
    //const uint MAX_L = maxLength(masks);
    LmerHeap lmers;

    cout <<"Scanning for significant Lmers..." << endl;

    // add any lmer that meets minimum for freq to lmers
    setLmers(sequence, seqLength, mask, MIN, L, lmers);
    //vector<LmerVector*> toReinsert;
    //for(int i = 0; i < 10; i++){
    //    cout << lmers.top()->size() << "\t" << lmers.top()->front() << endl;
    //    LmerVector* u = lmers.top();
    //    lmers.pop();
    //    toReinsert.push_back(u);
    //}
    //for(vector<LmerVector*>::iterator it = toReinsert.begin(); it != toReinsert.end(); ++it) {
    //    lmers.push(*it);
    //}

    cout <<"Finding elementary repeat families..." <<endl;

    // go through LmerHeap one at a time (sorted by decreasing freq/increasing position)
    while(!lmers.empty()){
        createNewFamily(lmers, L, families);
    }
    //cout <<"Splitting incomplete families..." <<endl;
}

#endif //SCANER_SEEDCHAIN
