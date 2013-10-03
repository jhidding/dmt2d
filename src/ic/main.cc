#include "../base/system.hh"
#include "../base/format.hh"
#include "ic.hh"

#include <fstream>

using namespace System;

void cmd_ic(int argc, char **argv)
{
	std::ostringstream ss;
	ss << time(NULL);
	std::string timed_seed = ss.str();

	Argv C = read_arguments(argc, argv,
		Option({0, "h", "help", "false",
			"print help on the use of this program."}),

		Option({Option::VALUED | Option::CHECK, "i", "id", date_string(),
			"identifier for filenames."}),

		Option({Option::VALUED | Option::CHECK, "d", "dim", "2",
			"number of dimensions in which to run the code. "
			"Usually this should by 2 or 3."}),

		Option({Option::VALUED | Option::CHECK, "b", "mbits", "8",
			"size of the mass-box, or binary power of the number "
		       	"of particles per dimension. A value of 8 (default) "
			"will give a boxsize of 256."}),

		Option({Option::VALUED | Option::CHECK, "L", "size", "100",
			"physical size of the simulation in units of Mpc/h."}),

		Option({Option::VALUED | Option::CHECK, "", "seed", timed_seed,
			"random seed used to generate the initial conditions. "
			"By default the number of seconds since Epoch is used. "}),

		Option({Option::VALUED | Option::CHECK, "n", "power-slope", "-1.0",
			"the slope of the power spectrum."}),

		Option({0, "", "smooth", "false",
			"whether to smooth the initial conditions with a Gaussian "
			"filter or not."}),

		Option({Option::VALUED | Option::CHECK, "", "scale", "4.0",
			"scale at which to smooth in units of Mpc/h."}),
		
		Option(0, "p", "potential", "false",
			"include the potential in the result."),
		
		Option(0, "z", "displacement", "false",
			"include the Zel'dovich displacement in the result." ));

	if (C.get<bool>("help"))
	{
		std::cout << "Cosmic workset Conan, by Johan Hidding.\n\n";
		C.print(std::cout);
		exit(0);
	}

	Header H; H << C;
	History I; I << C;
	Array<double> D = Conan::generate_random_field(H);

	std::ofstream fo(C["id"] + ".density.init.conan"); H.to_file(fo); I.to_file(fo);
	save_to_file(fo, D, "density");

	if (C.get<bool>("potential"))
	{
		Conan::compute_potential(H, D);
		save_to_file(fo, D, "potential");
	}

	if (C.get<bool>("displacement"))
	{
		Conan::compute_displacement(H, D, fo);
	}

	fo.close();
}

Global<Command> _IC("ic", cmd_ic);

