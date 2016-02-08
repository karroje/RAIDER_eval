#ifndef LMER_VECTOR_WRAPPER
#define LMER_VECTOR_WRAPPER
#include "LmerVector.h"

using namespace std;

#ifndef uint
typedef unsigned int uint;
#endif


class LmerVectorWrapper {
public:
    //LmerVectorWrapper(uint seed, uint count) {
    //    lmer = new LmerVector();
    //    this.seed = seed;
    //    this.count = count;
    //}

    LmerVectorWrapper(uint seed, vector<uint> positions) {
        lmer = new LmerVector();
        for(std::vector<T>::iterator it = positions.begin(); it != positions.end(); ++it) {
            lmer.push_back(*it);
        }
        this.seed = seed;
        this.count = positions.size();
        this.first = positions.front();
    }

    LmerVector* getLmer(){
        return lmer;
    }

    uint size(){
        return lmer.size();
    }


    uint seed;
    uint count;
    uint first;
    LmerVector* lmer;
};

#endif //LMER_VECTOR
