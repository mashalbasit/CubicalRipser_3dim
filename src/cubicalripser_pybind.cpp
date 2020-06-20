/*
This file is part of CubicalRipser
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
Modified by Shizuo Kaji

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <fstream>
#include <iostream>
#include <algorithm>
#include <queue>
#include <vector>
#include <unordered_map>
#include <string>
#include <cstdint>


#include "cube.h"
#include "dense_cubical_grids.h"
#include "coboundary_enumerator.h"
#include "write_pairs.h"
#include "joint_pairs.h"
#include "compute_pairs.h"
#include "config.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

using namespace std;

namespace py = pybind11;

/////////////////////////////////////////////
py::array_t<double> computePH(py::array_t<double> img, int maxdim=0, const std::string &location="birth"){

	Config config;
	config.format = NUMPY;

	vector<WritePairs> writepairs; // (dim birth death x y z)
	writepairs.clear();
	
	DenseCubicalGrids* dcg = new DenseCubicalGrids(config);
	vector<Cube> ctr;

    const auto &buff_info = img.request();
    const auto &shape = buff_info.shape;
	dcg->dim = buff_info.ndim;
	config.maxdim = maxdim;
	config.maxdim = std::min<uint8_t>(config.maxdim, dcg->dim - 1);
	if(location=="death") config.location=LOC_DEATH;

	// load file
	dcg->ax = shape[0];
	if (dcg->dim == 3) {
		dcg->ay = shape[1];
		dcg->az = shape[2];
		for (uint32_t x = 0; x < dcg->ax; ++x) {
			for (uint32_t y = 0; y < dcg->ay; ++y) {
				for (uint32_t z = 0; z < dcg->az; ++z) {
					dcg->setBirth(x,y,z,*img.data(x, y, z));
				}
			}
		}
	}else if (dcg->dim == 2) {
		dcg->ay = shape[1];
		dcg->az = 1;
		for (uint32_t x = 0; x < dcg->ax; ++x) {
			for (uint32_t y = 0; y < dcg->ay; ++y) {
				for (uint32_t z = 0; z < dcg->az; ++z) {
					dcg->setBirth(x,y,z,*img.data(x, y));
				}
			}
		}
	}else if (dcg->dim == 1) {
		dcg->ay = 1;
		dcg->az = 1;
		for (uint32_t x = 0; x < dcg->ax; ++x) {
			for (uint32_t y = 0; y < dcg->ay; ++y) {
				for (uint32_t z = 0; z < dcg->az; ++z) {
					dcg->setBirth(x,y,z,*img.data(x));
				}
			}
		}
	}

	dcg -> axy = dcg->ax * dcg->ay;
	dcg -> ayz = dcg->ay * dcg->az;
	dcg -> axyz = dcg->ax * dcg->ay * dcg->az;



	// compute PH
	ComputePairs* cp = new ComputePairs(dcg, writepairs, config);
    vector<uint32_t> betti(0);

	JointPairs* jp = new JointPairs(dcg, writepairs, config);
	if(dcg->dim==1){
		jp -> enum_edges({0},ctr);
	}else if(dcg->dim==2){
		jp -> enum_edges({0,1},ctr);
	}else{
		jp -> enum_edges({0,1,2},ctr);
	}
	jp -> joint_pairs_main(ctr,0); // dim0
	betti.push_back(writepairs.size());
	if(config.maxdim>0){
		cp -> compute_pairs_main(ctr); // dim1
		betti.push_back(writepairs.size() - betti[0]);
		if(config.maxdim>1){
			cp -> assemble_columns_to_reduce(ctr,2);
			cp -> compute_pairs_main(ctr); // dim2
			betti.push_back(writepairs.size() - betti[0] - betti[1]);
		}
	}

	// result
	int64_t p = writepairs.size();
	vector<ssize_t> result_shape{p,6};
	py::array_t<double> data{result_shape};
	for(int64_t i = 0; i < p; ++i){
		*data.mutable_data(i, 0) = writepairs[i].dim;
		*data.mutable_data(i, 1) = writepairs[i].birth;
		*data.mutable_data(i, 2) = writepairs[i].death;
		*data.mutable_data(i, 3) = writepairs[i].birth_x;
		*data.mutable_data(i, 4) = writepairs[i].birth_y;
		*data.mutable_data(i, 5) = writepairs[i].birth_z;
	}
	return data;
}

PYBIND11_MODULE(cripser, m) {
    m.doc() = R"pbdoc(
        Cubical Ripser plugin
        -----------------------
        .. currentmodule:: cripser
        .. autosummary::
           :toctree: _generate
           add
           subtract
    )pbdoc";

    m.def("computePH", &computePH, R"pbdoc(Compute Persistent Homology
    )pbdoc", py::arg("arr"),  py::arg("maxdim")=2, py::arg("location")="birth");

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}