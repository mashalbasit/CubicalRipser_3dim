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

using namespace std;

enum file_format { DIPHA, PERSEUS, NUMPY };

class DenseCubicalGrids {
public:
	double threshold;
	int dim;
	int ax, ay, az;
    long axy, axyz, ayz;
	double*** dense3;

	DenseCubicalGrids(const std::string& filename, double _threshold, file_format format);
	~DenseCubicalGrids();
	
	double getBirthday(int x, int y, int z, int cm, int dim);
    double getBirthday(vector<int> &loc, int dim);
    double getBirthday(long ind, int dim);

	// unique id for each simplex (unique only within a single dimension)
	long getIndex(int x, int y, int z, int m=0){
		return(x + y * ax + z * axy + m * axyz);
	}
    
    // conversion from id to coordinates and type
    vector<int> getXYZM(long index) {
        vector<int> loc(4);   // (x,y,z,m)
        loc[0] = index % ax;
        loc[1] = (index / ax) % ay;
        loc[2] = (index / axy) % az;
        loc[3] = (index / axyz);
        return(loc);
    }

};
