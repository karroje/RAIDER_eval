// ==========================================================================
//                               	RAIDER2
// ==========================================================================
// Author: Charlotte Schaeffer <schaefce@miamioh.edu>
// Based off of original RAIDER by Nathan Figueroa
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <unordered_map>

#include "SeedChain.h"
#include "Family.h"
#include "LmerVector.h"

using namespace std;

typedef pair<uint, seqan::CharString> idThreshold;

const string OUTPUT_ELEMENTARY_REPEATS_FILENAME = "elements";
const string OUTPUT_REPEAT_FAMILIES_FILENAME = "families";
const string OUTPUT_COMPOSITES_FILENAME = "composites";
const string OUTPUT_SUMMARY_FILENAME = "summary_info";

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
struct AppOptions {
	// Age.  0 -- curr, 1 -- old, 2 -- oldest, 
	int age;
	// Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- carly verbose
	int verbosity;
	// Minimum repeat length
	uint min;
	// Minimum number of repeats to be significant
	uint count;

	seqan::CharString seed;
	seqan::CharString sequence_file;
	seqan::CharString output_directory;

	AppOptions() :
			verbosity(1) {}
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------
seqan::ArgumentParser::ParseResult parseCommandLine(AppOptions & options, int argc, char const ** argv) {
	// Setup ArgumentParser.
	seqan::ArgumentParser parser("RAIDER2");
	// Set short description, version, and date.
	setShortDescription(parser, "RAIDER - Rapid Ab Initio Detection of Elementary Repeats");
	setVersion(parser, "2.0");
	setDate(parser, "June 2015");

	// Define usage line and long description.
	addUsageLine(parser, "[\\fIOPTIONS\\fP]  \"\\fISEQUENCE_FILE\\fP\"  \"\\fIOUTPUT_DIRECTORY\\fP\"");
	addDescription(
			parser,
			"RAIDER2 parses the given sequence file using the supplied mask (spaced seed) to identify de novo repeats. Minimum repeat size and other options can be configured as described below.");

	// We require two arguments.
	//addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "MASK"));
	addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "SEQUENCE_FILE"));
	addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "OUTPUT_DIRECTORY"));
	
    addOption(
			parser,
			seqan::ArgParseOption("s", "seed", "Spaced seed/mask to use. Defaults to 111110011111110001111111000000000000011111.",
					seqan::ArgParseOption::STRING));

	addOption(
			parser,
			seqan::ArgParseOption("m", "min", "Minimum repeat length. Defaults to pattern length.",
					seqan::ArgParseOption::INTEGER));
	addOption(
			parser,
			seqan::ArgParseOption("c", "count", "Minimum number of repeats in a family. Defaults to 5.",
					seqan::ArgParseOption::INTEGER));
	addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
	addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
	addOption(parser, seqan::ArgParseOption("a", "verbose+", "Enable extremely verbose output."));
	addOption(
			parser,
			seqan::ArgParseOption("g", "age", "Age of raiderv2. Defaults to 0.",
					seqan::ArgParseOption::INTEGER));

	// Add Examples Section.
	addTextSection(parser, "Examples");
	addListItem(parser, "\\fBraider\\fP \\fB-v\\fB -s \\fI1110110111\\fP \\fIchr23.fasta\\fP \"\\fIchr23_out\\fP\"",
			"Call with mask \"1110110111\" and verbose output.");

	// Parse command line.
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// Only extract  options if the program will continue after parseCommandLine()
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	// Extract option values.
	if (isSet(parser, "quiet"))
		options.verbosity = 0;
	if (isSet(parser, "verbose"))
		options.verbosity = 2;
    if (isSet(parser, "verbose+"))
        options.verbosity = 3;

	//seqan::getArgumentValue(options.seed, parser, 0);
	seqan::getArgumentValue(options.sequence_file, parser, 0);
	seqan::getArgumentValue(options.output_directory, parser, 1);


	if (isSet(parser, "count"))
		seqan::getOptionValue(options.count, parser, "count");
	else
		options.count = 5;

    if (isSet(parser, "seed"))
        seqan::getOptionValue(options.seed, parser, "seed");
    else
        options.seed = "111110011111110001111111000000000000011111";

	if (isSet(parser, "min"))
		seqan::getOptionValue(options.min, parser, "min");
	else
		options.min = seqan::length(options.seed);
	
    if (isSet(parser, "age"))
		seqan::getOptionValue(options.age, parser, "age");
	else
		options.age = 0;
	// Ensure a trailing /
	if (options.output_directory[seqan::length(options.output_directory) - 1] != '/') {
		seqan::append(options.output_directory, "/");
	}

	return seqan::ArgumentParser::PARSE_OK;
}

/**
 * Given a sequence stream, read in all sequences and concatenate into a master sequence
 */
bool concatenateSequences(seqan::SequenceStream &seqStream, vector<idThreshold> &thresholds,
		seqan::Dna5String &outSequence, int verbosity) {
	seqan::CharString id;
	if (seqan::readRecord(id, outSequence, seqStream) != 0) {
		return false;
	}
	if (verbosity > 0) {
		cout << "Preparing " << id << endl;
	}
	uint currentThreshold = seqan::length(outSequence);
	thresholds.push_back(make_pair(currentThreshold, id));
	seqan::Dna5String other;
	while (seqan::readRecord(id, other, seqStream) == 0) {
		if (verbosity > 0) {
			cout << "Preparing " << id << endl;
		}
		currentThreshold += seqan::length(other);
		thresholds.push_back(make_pair(currentThreshold, id));
		seqan::append(outSequence, other);
	}
	return true;
}

/**
 * Given a file name, read all sequences from file and compile into master sequence. Also return
 * list of thresholds: indices into the master sequence paired with sequence IDs to mark where each
 * individual sequence begins.
 */
bool getSequence(seqan::CharString file, seqan::Dna5String &sequence, vector<idThreshold> &thresholds, int verbosity) {
	if (verbosity > 0) {
		cout << "Loading sequence..." << endl;
	}

	seqan::SequenceStream seqStream(seqan::toCString(file));
	if (!seqan::isGood(seqStream) || !concatenateSequences(seqStream, thresholds, sequence, verbosity)) {
		cout << "Error: unable to open sequence." << endl;
		return false;
	}
	return true;
}

/**
 * Given an index into the master sequence, return the associated sequence ID
 */
seqan::CharString getSeqId(uint index, vector<idThreshold> &thresholds) {
	for (uint i = 0; i < thresholds.size(); i++) {
		if (index < thresholds[i].first) {
			return thresholds[i].second;
		}
	}
	cout << "ERROR: unable to find associated sequence ID" << endl;
	return "Unknown";
}

/**
 * Given an index into the master sequence, find the corresponding original sequence and
 * return the corresponding index into that sequence
 */
uint getSeqIndex(uint index, vector<idThreshold> &thresholds) {
	uint prevThreshold = 0;
	for (uint i = 0; i < thresholds.size(); i++) {
		if (index < thresholds[i].first) {
			return index - prevThreshold;
		}
		prevThreshold = thresholds[i].first;
	}
	cout << "ERROR: unable to find correct threshold" << endl;
	// TODO throw exceptions
	return index;
}

/**
 * Print the command line arguments back to the user.
 */
void printArgs(AppOptions &options) {
	if (options.verbosity > 0) {
		cout << "__ARGUMENTS____________________________________________________________________" << endl
				<< "VERBOSITY\t" << options.verbosity << endl
				<< "MIN_LENGTH\t" << options.min << endl
				<< "MIN_COUNT\t" << options.count << endl
				<< "SPACED_SEED     \t" << options.seed <<endl
				<< "SEQUENCE_FILE\t" << options.sequence_file << endl
				<< "OUTPUT_DIRECTORY\t" << options.output_directory << endl;
	}
}

void writeFamilies(vector<Family*> &families, AppOptions &options) {
	ofstream ofile;
	string familiesFilePath(seqan::toCString(options.output_directory));
	familiesFilePath.append(OUTPUT_REPEAT_FAMILIES_FILENAME);

	if (options.verbosity > 0) {
		cout << "Writing families to " << familiesFilePath << endl;
	}

	ofile.open(familiesFilePath.c_str());

	ofile << "#  fam\tnum_copies\tcopy_length" << endl;
	for (uint i = 0; i < families.size(); i++) {
		Family* fam = families[i];
		if (fam->repeatLength(options.min) >= options.min && fam->size() >= options.count) {
			ofile << i << "\t" << fam->size() << "\t" << fam->repeatLength(options.min) << endl;
		}
	}
}

void writeRepeats(vector<Family*> &families, vector<idThreshold> &thresholds, AppOptions &options) {
	uint repCount = 0;

	ofstream ofile;
	string repeatsFilePath(seqan::toCString(options.output_directory));
	repeatsFilePath.append(OUTPUT_ELEMENTARY_REPEATS_FILENAME);

	if (options.verbosity > 0) {
		cout << "Writing elementary repeats to " << repeatsFilePath << endl;
	}

	ofile.open(repeatsFilePath.c_str());

	ofile << "#  fam\tele\tdir\tsequence\tstart\tend" << endl;

	for (uint i = 0; i < families.size(); i++) {
		Family* fam = families[i];

		if (fam->repeatLength(options.min) >= options.min && fam->size() >= options.count) {
			int famId = i;
			int repId = repCount;
			repCount += fam->size();

			LmerVector* prefix = fam->getPrefix();
			for (uint i = 0; i < prefix->size(); i++) {
				repId++;
				uint index = (*prefix)[i];
				uint trueIndex = getSeqIndex(index, thresholds);
				uint length = fam->repeatLength(options.min);
				ofile << famId << "\t" << repId << "\t" << "1" << "\t" << getSeqId(index, thresholds) << "\t"
						<< trueIndex << "\t" << trueIndex + length << endl;
			}
		}

	}
}


void writeSummary(vector<Family*> &families, AppOptions &options) {
	uint maxLen = 0;
	uint maxSize = 0;
	uint maxSecondSize = 0;
	uint repCount = 0;
	uint famCount = 0;

	ofstream ofile;
	string summaryFilePath(seqan::toCString(options.output_directory));
	summaryFilePath.append(OUTPUT_SUMMARY_FILENAME);

	if (options.verbosity > 0) {
		cout << "Writing summary to " << summaryFilePath << endl;
	}

	ofile.open(summaryFilePath.c_str());

	for (uint i = 0; i < families.size(); i++) {
		Family* fam = families[i];
		if (fam->repeatLength(options.min) >= options.min && fam->size() >= options.count) {
			famCount++;
			repCount += fam->size();

			if (fam->size() > maxSize) {
				maxSecondSize = maxSize;
				maxSize = fam->size();
			}

			if (fam->repeatLength(options.min) > maxLen) {
				maxLen = fam->repeatLength(options.min);
			}
		}
	}

	ofile << "#FAMILIES\t" << famCount << endl;
	ofile << "#REPEATS\t" << repCount << endl;
	ofile << "LONGEST_REPEAT\t" << maxLen << endl;
	ofile << "LONGEST_FAMILY\t" << maxSize << endl;
	ofile << "SECOND_LONGEST\t" << maxSecondSize << endl;
}

void writeResults(vector<Family*> &families, vector<idThreshold> &thresholds, AppOptions &options) {
	writeFamilies(families, options);
	writeRepeats(families, thresholds, options);
	writeSummary(families, options);
}


int main(int argc, char const ** argv) {
	AppOptions options;
	parseCommandLine(options, argc, argv);
	printArgs(options);

	// Load sequences into master sequence, mark where each sequence ID begins within the master
	vector<idThreshold> thresholds;
	seqan::Dna5String sequence;
	if (getSequence(options.sequence_file, sequence, thresholds, options.verbosity) == false) {
		// TODO throw exceptions
		return 1;
	}

	if (options.verbosity > 0) {
		cout << "BASE PAIRS\t" << seqan::length(sequence) << endl;
	}

	vector<Family*> families;
	vector<seqan::CharString> seeds;
	seeds.push_back(options.seed);

	getElementaryFamilies(sequence, seeds, families, options.verbosity, options.age);

	if (options.verbosity > 0) {
		cout << "Writing results elements..." << endl;
	}
	writeResults(families, thresholds, options);

	return 0;
}
