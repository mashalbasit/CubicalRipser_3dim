/* dense_cubical_grids.cpp

This file is part of CubicalRipser
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
Modified by Shizuo Kaji

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cassert>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
#include <initializer_list>

#include "dense_cubical_grids.h"
#include "npy.hpp"

using namespace std;


DenseCubicalGrids::DenseCubicalGrids(Config& _config)  {
	config = &_config;
	threshold = config->threshold;
}


// read from file
void DenseCubicalGrids::loadImage(bool _embedded){
	embedded = _embedded;
	// read file
	cout << "Reading " << config->filename << endl;
	switch(config->format){
		case DIPHA:
		{
			ifstream fin( config->filename, ios::in | ios::binary );

			int64_t d;
			fin.read( ( char * ) &d, sizeof( int64_t ) ); // magic number
			assert(d == 8067171840);
			fin.read( ( char * ) &d, sizeof( int64_t ) ); // type number
			assert(d == 1);
			fin.read( ( char * ) &d, sizeof( int64_t ) ); //data num
			fin.read( ( char * ) &d, sizeof( int64_t ) ); // dim 
			dim = d;
			assert(dim < 4);
			fin.read( ( char * ) &d, sizeof( int64_t ) );
			ax = d;
			if (dim>1) {
				fin.read( ( char * ) &d, sizeof( int64_t ) );
				ay = d;
			}else{
				ay = 1;
			}
			if (dim>2) {
				fin.read((char *)&d, sizeof(int64_t));
				az = d;
			}else {
				az = 1;
			}
			cout << "ax : ay : az = " << ax << " : " << ay << " : " << az << endl;
            data.resize(ax*ay*az);

			double dou;
			for (uint32_t z = 0; z < az; ++z) {
				for (uint32_t y = 0; y < ay; ++y) {
					for (uint32_t x = 0; x < ax; ++x) {
                        if (!fin.eof()) {
                            fin.read((char *)&dou, sizeof(double));
                            setBirth(x,y,z,dou);
                        }
                        else {
                            cerr << "file endof error " << endl;
                        }
					}
				}
			}
			fin.close();
			break;
		}

		case PERSEUS:
		{
			ifstream reading_file; 
			reading_file.open(config->filename.c_str(), ios::in); 

			string reading_line_buffer; 
			getline(reading_file, reading_line_buffer); 
			dim = atoi(reading_line_buffer.c_str());
			assert(dim < 4);
			getline(reading_file, reading_line_buffer);
			ax = atoi(reading_line_buffer.c_str()); 
			if (dim>1) {
				getline(reading_file, reading_line_buffer); 
				ay = atoi(reading_line_buffer.c_str()); 
			}else {
				ay = 1;
			}
			if (dim>2) {
				getline(reading_file, reading_line_buffer);
				az = atoi(reading_line_buffer.c_str());
			}else {
				az = 1;
			}
			cout << "ax : ay : az = " << ax << " : " << ay << " : " << az << endl;
            data.resize(ax*ay*az);
			for (uint32_t z = 0; z < az; ++z) {
				for (uint32_t y = 0; y < ay; ++y) {
					for (uint32_t x = 0; x < ax; ++x) {
						if (!reading_file.eof()) {
							getline(reading_file, reading_line_buffer);
							double dou = atof(reading_line_buffer.c_str());
							if (dou == -1) {
								setBirth(x,y,z,config->threshold);
							}
							else {
								setBirth(x,y,z,dou);
							}
						}
					} 
				}
			}
			reading_file.close();
			break;
		}
		case NUMPY:
		{
			vector<unsigned long> shape;
			npy::LoadArrayFromNumpy(config->filename.c_str(), shape, data);
			if(shape.size() > 3){
				cerr << "Input array should be 2 or 3 dimensional " << endl;
				exit(-1);
			}
			dim = shape.size();
			ax = shape[0];
			if (dim>1) {
				ay = shape[1];
			}else {
				ay = 1;
			}
			if (dim>2) {
				az = shape[2];
			}else {
				az = 1;
			}
			cout << "ax : ay : az = " << ax << " : " << ay << " : " << az << endl;
			break;
		}
	}	
    if(embedded){ // dual complex
        if(az>1){
            az = az + 2;
        }
        ax = ax + 2;
        ay = ay + 2;
    }
	axy = ax * ay;
	ayz = ay * az;
	axyz = ax * ay * az;
//	cout << ax << "," << ay << "," << az << endl;
}

// return filtlation value for a cube

double DenseCubicalGrids::getBirth(int32_t x, int32_t y, int32_t z, uint8_t cm, uint8_t dim) {
	// beware of the shift due to the boundary
	switch (dim) {
		case 0:
            return(getBirth(x,y,z));
		case 1:
			switch (cm) {
			case 0:
				return max(getBirth(x,y,z), getBirth(x+1,y  ,z));
			case 1:
				return max(getBirth(x,y,z), getBirth(x  ,y+1,z));
			case 2:
				return max(getBirth(x,y,z), getBirth(x  ,y  ,z+1));
			case 3:
				return max(getBirth(x,y,z), getBirth(x+1,y+1,z));
			case 4:
				return max(getBirth(x,y,z), getBirth(x+1,y-1,z));
			// for 3d dual only
			case 5:
				return max(getBirth(x,y,z), getBirth(x  ,y-1,z+1));
			case 6:
				return max(getBirth(x,y,z), getBirth(x  ,y+1,z+1));
			case 7:
				return max(getBirth(x,y,z), getBirth(x+1,y-1,z+1));
			case 8:
				return max(getBirth(x,y,z), getBirth(x+1,y  ,z+1));
			case 9:
				return max(getBirth(x,y,z), getBirth(x+1,y+1,z+1));
			case 10:
				return max(getBirth(x,y,z), getBirth(x+1,y-1,z-1));
			case 11:
				return max(getBirth(x,y,z), getBirth(x+1,y  ,z-1));
			case 12:
				return max(getBirth(x,y,z), getBirth(x+1,y+1,z-1));
			}
		case 2:
			switch (cm) {
			case 0: // x - y (fix z)
				return max({ getBirth(x,y,z), getBirth(x+1,y,z),
					getBirth(x+1,y+1,z), getBirth(x,y+1,z) });
			case 1: // z - x (fix y)
				return max({ getBirth(x,y,z), getBirth(x,y,z+1),
					getBirth(x+1,y,z+1), getBirth(x+1,y,z) });
			case 2: // y - z (fix x)
				return max({ getBirth(x,y,z), getBirth(x,y+1,z),
					getBirth(x,y+1,z+1), getBirth(x,y,z+1) });
			}
		case 3:
			return max({ getBirth(x,y,z), getBirth(x+1,y,z),
				getBirth(x+1,y+1,z), getBirth(x,y+1,z),
				getBirth(x,y,z+1), getBirth(x+1,y,z+1),
				getBirth(x+1,y+1,z+1), getBirth(x,y+1,z+1) });
		}
	return threshold; // dim > 3
}

// allocate 3d array
double ***DenseCubicalGrids::alloc3d(uint32_t x, uint32_t y, uint32_t z) {
	double ***d = (double***)malloc(x * sizeof(double**));
	d[0] = (double**)malloc(x * y * sizeof(double*));
	d[0][0] = (double*)malloc(x*y*z * sizeof(double));
	for (uint32_t i = 0; i < x ; i++) {
		d[i] = d[0] + i * y;
		for (uint32_t j = 0; j < y; j++) d[i][j] = d[0][0] + i * y*z + j * z;
	}
	if (d == NULL) {
		cerr << "not enough memory!" << endl;
	}
	return d;
}

DenseCubicalGrids::~DenseCubicalGrids(){
//	free(dense3[0][0]);
//	free(dense3[0]);
//	free(dense3);
}
