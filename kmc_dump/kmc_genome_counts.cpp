#include "stdafx.h"
#include <iostream>
#include <fstream>
#include "./../kmc_api/kmer_api.h"
#include "./../kmc_api/kmc_file.h"
#include "./../kmc_tools/defs.h"
#include "nc_utils.h"


class ChromInfo
{
public:
	int kmer_size;
	uint32_t max_counter;
	std::string id;
	int size;
	std::vector<std::vector<uint32_t>> *scores;
	int chunks;
	char *buffer;
	int start;
	int end;
	bool last_chrom;
	int N;
	int N0;
	ChromInfo(int _size, uint32_t _max_counter);
	void ExpandChunks();
	void AddChunk();
	void RotateBuffer();
	void ScoreKmers(CKMCFile *db);
	void WriteWiggle();
};

class ChromFastaReader
{
	bool end_reached;
    std::fstream in;
    bool opened;
    std::string single;
public:
	ChromFastaReader();
	~ChromFastaReader();
	void SetSingle(std::string _single);
	bool OpenFile(std::string _input_file_name);
	int GetChunk(ChromInfo *chrom);
	bool Finished(); 
};

void print_info(void);

//----------------------------------------------------------------------------------
// Check if --help or --version was used
bool help_or_version(int argc, char** argv)
{
	const std::string version = "--version";
	const std::string help = "--help";
	for (int i = 1; i < argc; ++i)
	{
		if (argv[i] == version || argv[i] == help)
			return true;
	}
	return false;
}

int _tmain(int argc, char* argv[])
{
	if (argc == 1 || help_or_version(argc, argv))
	{
		print_info();
		return 0;
	}

	CKMCFile *kmer_data_base;
	int32 i;
	uint32_t max_counter=255;
	std::string single="";
	std::string input_file_name;
	std::string fasta_file_name;

	//------------------------------------------------------------
	// Parse input parameters
	//------------------------------------------------------------
	if(argc < 3)
	{
		print_info();
		return 0;
	}

	for(i = 1; i < argc; ++i)
		if(argv[i][0] == '-')
		{
			if(strncmp(argv[i], "-cx", 3) == 0)
				max_counter = atoi(&argv[i][3]);
			else if (strncmp(argv[i], "-ch", 3) == 0)
				single.assign(&argv[i][3]);
		}
		else
			break;

	if(argc - i < 2)
	{ 
		print_info();
		return 0;
	}

	input_file_name = std::string(argv[i++]);
	fasta_file_name = std::string(argv[i]);

	kmer_data_base = new CKMCFile;
	kmer_data_base->OpenForRA(input_file_name);
	int kmer_len = kmer_data_base->KmerLength();

	ChromFastaReader reader;
	reader.SetSingle(single);
	reader.OpenFile(fasta_file_name);
	ChromInfo chrom(kmer_len, max_counter);
	while (true){
		chrom = ChromInfo(kmer_len, max_counter);
		bool finished=false;
		while (!finished)
		{
			finished = reader.GetChunk(&chrom);
			if (chrom.id.length() == 0)
				return 0;
			std::cerr << "\r" << chrom.id << ":" << chrom.N0 * chrom.chunks;
			chrom.ScoreKmers(kmer_data_base);
			if (!finished)
				chrom.RotateBuffer();
		}
		std::cerr << "\r                            \r";
		std::cerr << chrom.id << " (" << chrom.size << ")\n";
		chrom.WriteWiggle();
		if (chrom.last_chrom)
			break;
	}
	return 1;
}

// -------------------------------------------------------------------------
// Print execution options 
// -------------------------------------------------------------------------
void print_info(void)
{
	std::cout << "KMC get genome counts ver. " << KMC_VER << " (" << KMC_DATE << ")\n"
			  << "\nUsage:\nget_genome_counts [options] <kmc_database> <genome_fasta_file>\n"
			  << "Parameters:\n"
			  << "<kmc_database>      - kmer_counter's output\n"
			  << "<genome_fasta_file> - fasta file used to generate kmc_database\n"
			  << "Options:\n"
			  << "-cx<value> - truncate k-mer counts occurring more of than <value> times (default 255)\n"
			  << "-ch<value> - score only chromosome <value>\n";
}

/************************************************
ChromInfo Class Functions
************************************************/

ChromInfo::ChromInfo(int _size, uint32_t _max_counter)
{
	N0 = 1000000;
	N = N0 + _size - 1;
	kmer_size = _size;
	max_counter = _max_counter;
	buffer = new char[N];
	id = "";
	size = 0;
	chunks = 0;
	scores = new std::vector<std::vector<uint32_t>>;
	start = 0;
	last_chrom = false;
}

void ChromInfo::AddChunk()
{
	std::vector<uint32_t> new_scores;
	std::vector<bool> new_valid;
	scores->push_back(new_scores);
	chunks += 1;
}

void ChromInfo::RotateBuffer()
{
	for (int i=0; i<kmer_size-1; i++)
		buffer[i] = buffer[N0 + i];
	start = kmer_size - 1;
}

void ChromInfo::ScoreKmers(CKMCFile *db)
{
	std::string seq="";
	for (int i=0; i<end; i++)
		seq += buffer[i];
	std::vector<uint32_t> *tmp_scores = &(scores->at(chunks-1));
	db->GetCountersForRead(seq, *tmp_scores);

	int n_count=0;
	char cN='N';
	for (int i=0; i<kmer_size; i++)
		if (buffer[i] == cN)
			n_count ++;
	if (n_count > 0)
		tmp_scores->at(0) = 0;
	else
		if (tmp_scores->at(0) == 0)
			tmp_scores->at(0) = 1;
		else
			if (tmp_scores->at(0) > max_counter)
				tmp_scores->at(0) = max_counter;
	for (uint i=1; i<tmp_scores->size(); i++)
	{
		if (buffer[i-1] == cN)
			n_count --;
		if (buffer[i+kmer_size-1] == cN)
			n_count ++;
		if (n_count > 0)
			tmp_scores->at(i) = 0;
		else if (tmp_scores->at(i) == 0)
			tmp_scores->at(i) = 1;
		else if (tmp_scores->at(i) > max_counter)
			tmp_scores->at(i) = max_counter;
	}
}

void ChromInfo::WriteWiggle()
{
	std::string tmp_str="";
	int pos=1, chunk=0, i, stop;
	uint32_t score;
	while (chunk < chunks)
	{
		std::vector<uint32_t> *tmp_scores = &(scores->at(chunk));
		chunk ++;
		i = 0;
		stop = tmp_scores->size();
		while (i < stop)
		{
			score = tmp_scores->at(i);
			if (score > 0)
			{
				if (tmp_str.length() == 0)
				{
					tmp_str = "fixedStep chrom=" + id + " start=" + std::to_string(pos) + " step=1\n";
					//outfile->write(tmp_str.c_str(), tmp_str.length());
					std::cout << tmp_str;
				}
				tmp_str = std::to_string(score);
				//outfile->write(tmp_str.c_str(), tmp_str.length());
				//outfile->write("\n", 1);
				std::cout << tmp_str << "\n";
			}
			else
				tmp_str = "";
			pos += 1;
			i += 1;
		}
	}
}

/************************************************
ChromFastaReader Class Functions
************************************************/

// Constructor of FASTA reader
ChromFastaReader::ChromFastaReader()
{
	end_reached = false;
	opened = false;
	single = "";
}

// Destructor - close the files
ChromFastaReader::~ChromFastaReader()
{
	if(opened)
		in.close();
}

void ChromFastaReader::SetSingle(std::string _single)
{
	single.assign(_single);
}

// Open the file
bool ChromFastaReader::OpenFile(std::string _input_file_name)
{
	if(opened)
		return false;

	in.open(_input_file_name.c_str(), std::fstream::in);
	if(!in.is_open())
		return false;
	char test, gt='>';
	in.get(test);
	if (test != gt)
	{
		std::cerr << "Doesn't appear to be a fasta file\n";
		return false;
	}
	in.seekg(0);
	opened = true;

	return true;
}

int ChromFastaReader::GetChunk(ChromInfo *chrom)
{
	if (end_reached)
	{
		std::cerr << "Tried to read past last chromosome in FASTA file";
		return false;
	}

	chrom->AddChunk();

	char temp=' ', nl='\n', gt='>';
	long i=chrom->start, pos;
	bool finished=false;
	std::string id="";

	if (chrom->id.length() == 0)
	{
		while (true)
		{
			while (temp != gt)
				if(!in.get(temp))
					return true;
			if(!in.get(temp))
				return true;
			while (temp != nl)
			{
				id += temp;
				if(!in.get(temp))
					return true;
			}
			if ((single.length() == 0) | (single.compare(id) == 0))
				break;
			id = "";
		}
		chrom->id = id;
	}

	chrom->end = chrom->N;
	while (i < chrom->N)
	{
		if (in.get(temp))
		{
			if (temp == nl)
				continue;
			if (temp == gt)
			{
				finished = true;
				chrom->end = i;
				pos = in.tellg();
				in.seekg(pos - 1);
				if (single.length() != 0)
					chrom->last_chrom = true;
				break;
			}
			chrom->buffer[i] = toupper(temp);
		}
		else
		{
			finished = true;
			chrom->last_chrom = true;
			chrom->end = i;
			break;
		}
		i++;
	}
	chrom->size += chrom->end - chrom->start;
	return finished;
}

bool ChromFastaReader::Finished()
{
	if (end_reached)
		return true;
	else
		return false;
}
