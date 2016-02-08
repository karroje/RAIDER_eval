#include "gtest/gtest.h"
#include "Raider.h"
#include "Family.h"
#include <vector>
#include <stdexcept>


using namespace std;


TEST(LmerVector, setFamily_getFamily) {
	LmerVector r(0);
	EXPECT_EQ(r.getFamily(), nullptr);
	Family fam;
	fam.adopt(&r);
	r.setFamily(&fam);
	EXPECT_EQ(r.getFamily(), &fam);
}

TEST(LmerVector, access_assign_push_back) {
	LmerVector r(5);
	r.push_back(7);
	EXPECT_EQ(r[0], 5);
	EXPECT_EQ(r[1], 7);
	r[1] = 9;
	EXPECT_EQ(r[1], 9);
}

TEST(Family, size) {
	LmerVector v(0);
	v.push_back(7);

	Family fam;
	fam.adopt(&v);
	EXPECT_EQ(fam.size(), 2);
}


TEST(SeedChain, baseToInt) {
	EXPECT_NO_THROW({
		size_t result = baseToInt('c');
		EXPECT_EQ(result, 1);
		result = baseToInt('G');
		EXPECT_EQ(result, 2);
	});

	ASSERT_THROW(baseToInt('z'), std::out_of_range);
}


TEST(SeedChain, seedToInt) {
	const char* seed = "ACGTAC";
	const uint L = 6;
	size_t iSeed = seedToInt(seed, L);

	EXPECT_EQ(iSeed, 433);
}

TEST(SeedChain, getNextSeed) {
	seqan::Dna5String sequence = "NNNNNNACGTNCGGTANNNNN";

	const uint seqLength = seqan::length(sequence);
	const uint mLength = 4;
	uint index = 0;

	char seed[5] = {0};
	EXPECT_TRUE(getNextSeed(sequence, index, seqLength, mLength, seed));
	EXPECT_TRUE((seed[0] == 'A' && seed[1] == 'C' && seed[2] == 'G' && seed[3] == 'T'));

	EXPECT_TRUE(getNextSeed(sequence, ++index, seqLength, mLength, seed));
	EXPECT_TRUE((seed[0] == 'C' && seed[1] == 'G' && seed[2] == 'G' && seed[3] == 'T'));

	EXPECT_TRUE(getNextSeed(sequence, ++index, seqLength, mLength, seed));
	EXPECT_TRUE((seed[0] == 'G' && seed[1] == 'G' && seed[2] == 'T' && seed[3] == 'A'));

	EXPECT_FALSE(getNextSeed(sequence, ++index, seqLength, mLength, seed));
}

TEST(SeedChain, setLmers) {
	seqan::Dna5String sequence = "ACGTNNNNNNACGTNCCGTANNNNN";

	const uint seqLength = seqan::length(sequence);
	const uint L = 3;
	LmerMap lmers;

	const uint MIN = 3;

	size_t acg = seedToInt("ACG", L);
	size_t cgt = seedToInt("CGT", L);

	setLmers(sequence, seqLength, MIN, L, lmers);
	EXPECT_EQ(lmers.count(acg), 0);
	EXPECT_EQ(lmers.count(cgt), 1);
}

TEST(SeedChain, getAndInsert) {
	LmerMap lmers;
	LmerVector v;

	size_t Lmer = 7;

	lmers[Lmer] = &v;

	EXPECT_EQ(getAndInsert(lmers, Lmer, 10), &v);
	EXPECT_EQ(v.back(), 10);
	EXPECT_EQ(v.front(), 10);

	LmerVector *u = getAndInsert(lmers, 9, 100);
	EXPECT_EQ(u, nullptr);
}

TEST(SeedChain, isNewFamily) {
	LmerVector v(10);
	v.push_back(100);

	Family* fam = 0;

	EXPECT_TRUE(isNewFamily(fam, &v));

	fam = new Family;
	fam->adopt(&v);

	LmerVector u(11);
	u.push_back(101);

	EXPECT_FALSE(isNewFamily(fam, &u));
	delete fam;

	Family newFam;
	newFam.adopt(&v);
	v.push_back(101);

	EXPECT_TRUE(isNewFamily(&newFam, &u));
}

TEST(Family, getPrefix) {
	LmerVector v(0);
	v.push_back(7);

	Family fam;
	fam.adopt(&v);
	EXPECT_EQ(fam.getPrefix(), &v);
}

TEST(Family, adopt) {
	LmerVector v(0);
	v.push_back(7);

	Family fam;
	fam.adopt(&v);

	LmerVector u(3);
	u.push_back(9);

	fam.adopt(&u);

	EXPECT_EQ(u.getFamily(), v.getFamily());
	EXPECT_EQ(fam.getLmers()->back(), &u);
}



TEST(SeedChain, splitRepeatsByLmer_keepV) {
	LmerVector v(0);
	v.push_back(10);

	Family fam;
	fam.adopt(&v);

	LmerVector u(1);
	u.push_back(11);
	fam.adopt(&u);

	LmerVector w(2);
	w.push_back(12);
	fam.adopt(&w);

	LmerVector x(3);
	x.push_back(23);
	fam.adopt(&x);

	Family* newFam = splitRepeatsByLmer(&fam, &w, true);
	EXPECT_EQ(v.getFamily(), &fam);
	EXPECT_EQ(u.getFamily(), &fam);
	EXPECT_EQ(w.getFamily(), &fam);
	EXPECT_EQ(x.getFamily(), newFam);

	EXPECT_EQ(fam.getSuffix(), &w);
	delete newFam;
}



TEST(SeedChain, splitRepeatsByLmer_no_keepV) {
	LmerVector v(0);
	v.push_back(10);

	Family fam;
	fam.adopt(&v);

	LmerVector u(1);
	u.push_back(11);
	fam.adopt(&u);

	LmerVector w(2);
	w.push_back(22);
	fam.adopt(&w);

	LmerVector x(3);
	x.push_back(23);
	fam.adopt(&x);

	Family* newFam = splitRepeatsByLmer(&fam, &w, false);
	EXPECT_EQ(v.getFamily(), &fam);
	EXPECT_EQ(u.getFamily(), &fam);
	EXPECT_EQ(w.getFamily(), newFam);
	EXPECT_EQ(x.getFamily(), newFam);

	EXPECT_EQ(fam.getSuffix(), &u);
	delete newFam;
}


TEST(Family, repeatLength) {
	Family fam;

	LmerVector v(0);
	v.push_back(10);
	fam.adopt(&v);

	LmerVector u(1);
	u.push_back(11);
	fam.adopt(&u);

	EXPECT_EQ(fam.repeatLength(1), 2);
	EXPECT_EQ(fam.repeatLength(5), 6);
}


TEST(SeedChain, getElementaryFamilies1) {
	seqan::Dna5String sequence = "TAAACTAGGTCACTGTAAACTTGGTCACT";
	uint L = 3;

	vector<Family*> families;
	getElementaryFamilies(sequence, 2, L, families);

	ASSERT_EQ(families.size(), 3);

	//ACT
	Family *fam = families[0];
	ASSERT_EQ(fam->size(), 4);
	EXPECT_EQ(fam->getLmers()->size(), 1);
	EXPECT_EQ(fam->getLmers()->front()->front(), 3);

	//TAAAC
	fam = families[1];
	ASSERT_EQ(fam->size(), 2);
	EXPECT_EQ(fam->getLmers()->size(), 3);
	EXPECT_EQ(fam->getLmers()->front()->front(), 0);

	//GGTCAC
	fam = families[2];
	ASSERT_EQ(fam->size(), 2);
	EXPECT_EQ(fam->getLmers()->size(), 4);
	EXPECT_EQ(fam->getLmers()->front()->front(), 7);
}


TEST(SeedChain, getElementaryFamilies2) {
	seqan::Dna5String sequence = "CCACGTACTNACGTNCCACGTAANTACACGTA";
	uint L = 3;

	vector<Family*> families;
	getElementaryFamilies(sequence, 2, L, families);

	ASSERT_EQ(families.size(), 5);

	//ACTT
	Family *fam = families[0];
	ASSERT_EQ(fam->size(), 4);
	EXPECT_EQ(fam->getLmers()->size(), 2);
	EXPECT_EQ(fam->getLmers()->front()->front(), 2);

	//CCA
	fam = families[1];
	ASSERT_EQ(fam->size(), 2);
	EXPECT_EQ(fam->getLmers()->size(), 1);
	EXPECT_EQ(fam->getLmers()->front()->front(), 0);

	//GTA
	fam = families[2];
	ASSERT_EQ(fam->size(), 3);
	EXPECT_EQ(fam->getLmers()->size(), 1);
	EXPECT_EQ(fam->getLmers()->front()->front(), 4);

	//TAC
	fam = families[3];
	ASSERT_EQ(fam->size(), 2);
	EXPECT_EQ(fam->getLmers()->size(), 1);
	EXPECT_EQ(fam->getLmers()->front()->front(), 5);

	//CAC
	fam = families[4];
	ASSERT_EQ(fam->size(), 3);
	EXPECT_EQ(fam->getLmers()->size(), 1);
	EXPECT_EQ(fam->getLmers()->front()->front(), 1);
}


TEST(SeedChain, getElementaryFamilies3) {
	seqan::Dna5String sequence = "ACGTACCTGNTACCTGNACGTACNTACCTGNACGTACCTG";
	uint L = 3;

	vector<Family*> families;
	getElementaryFamilies(sequence, 2, L, families);

	EXPECT_EQ(families.size(), 3);

	//TAC
	Family *fam = families[0];
	ASSERT_EQ(fam->size(), 5);
	EXPECT_EQ(fam->getLmers()->size(), 1);
	EXPECT_EQ(fam->getLmers()->front()->front(), 3);

	//ACGTA
	fam = families[1];
	ASSERT_EQ(fam->size(), 3);
	EXPECT_EQ(fam->getLmers()->size(), 3);
	EXPECT_EQ(fam->getLmers()->front()->front(), 0);

	//ACCTG
	fam = families[2];
	ASSERT_EQ(fam->size(), 4);
	EXPECT_EQ(fam->getLmers()->size(), 3);
	EXPECT_EQ(fam->getLmers()->front()->front(), 4);
}
