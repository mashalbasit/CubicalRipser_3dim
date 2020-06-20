/* dense_cubical_grids.h

This file is part of CubicalRipser
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
Modified by Shizuo Kaji

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <string>
#include <vector>
#include "config.h"

using namespace std;

class DenseCubicalGrids {
public:
	Config *config;
	double threshold;
	uint8_t dim;
	uint32_t ax, ay, az;
    uint32_t axy, axyz, ayz;
	bool embedded; // alexander dual (sphere embedding)
    vector<double> data;
//	double*** dense3;

	DenseCubicalGrids(Config&);
	~DenseCubicalGrids();
	void loadImage(bool embedded);
	double ***alloc3d(uint32_t x, uint32_t y, uint32_t z);
	double getBirth(int32_t x, int32_t y, int32_t z, uint8_t cm, uint8_t dim);
    inline double getBirth(int32_t x, int32_t y, int32_t z) {
        if(embedded){
            if (1 <= x && x < ax-1 &&  1 <= y && y < ay-1 && 1 <= z && z < az-1) {
                return -data[(z-1) + (y-1)*(az-2) + (x-1)*(ay-2)*(az-2)];
            }else if(-1 == x || x == ax || -1 == y || y == ay || z==-1 || z==az) {
                return(threshold);
            }else{
                return(-threshold);
            }
        }else{
            if (0 <= x && x < ax &&  0 <= y && y < ay && 0 <= z && z < az) {
                return data[z+y*az+x*ayz];
            }else{
                return(threshold);
            }
        }
    }
    inline void setBirth(int32_t x, int32_t y, int32_t z, double val){
        data[z+y*az+x*ay*az] = val;
    }
};
