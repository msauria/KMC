#include "./../kmc_api/kmc_file.h"

void print_info(void)
{
	std::cout << "KMC kmers database kmer sum ver. " << KMC_VER << " (" << KMC_DATE << ")\n"
			  << "\nUsage:\nkmc_kmer_sum <kmc_database>\n"
			  << "Parameters:\n"
			  << "<kmc_database>      - kmer_counter's output\n";
}

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
		return 1;
	}

	std::string kmc_file_name;

	//------------------------------------------------------------
	// Parse input parameters
	//------------------------------------------------------------
	if(argc < 2)
	{
		print_info();
		return 1;
	}

	kmc_file_name = std::string(argv[1]);

	CKMCFile kmc_data_base;
	if (!kmc_data_base.OpenForListing(kmc_file_name))
	{
		std::cerr << "Unable to open " << kmc_file_name << " as KMC kmer database\n";
		return 1;
	}
	uint32 _kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count;
	uint64 _max_count, _total_kmers;
	kmc_data_base.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length,
				       _signature_len, _min_count, _max_count, _total_kmers);
	uint64 kmer_sum=0;
	CKmerAPI kmer_object(_kmer_length);
	if (_mode) // quake compatible mode
	{
		float counter;
		while (kmc_data_base.ReadNextKmer(kmer_object, counter))
			kmer_sum += uint32(counter);
	}
	else
	{
		uint32 counter;
		while (kmc_data_base.ReadNextKmer(kmer_object, counter))
			kmer_sum += counter;
	}
	kmc_data_base.Close();
	std::cout << kmer_sum << "\n";
	return 0;
}




