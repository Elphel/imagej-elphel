package com.elphel.imagej.tileprocessor;
/**
 **
 ** Conflict - Represent a "conflit" between connected supertiles
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  Conflict.java is free software: you can redistribute it and/or modify
 **  it under the terms of the GNU General Public License as published by
 **  the Free Software Foundation, either version 3 of the License, or
 **  (at your option) any later version.
 **
 **  This program is distributed in the hope that it will be useful,
 **  but WITHOUT ANY WARRANTY; without even the implied warranty of
 **  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **  GNU General Public License for more details.
 **
 **  You should have received a copy of the GNU General Public License
 **  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ** -----------------------------------------------------------------------------**
 **
 */

public class Conflict{
	int nsTile= -1;
	int start_layer;
	int end_layer;
	boolean [] start_dirs = new boolean[24];
	//  0 ...  7 - ortho-diagonal-ortho
	//  0 ...  3 from start_layer to  end_layer: 0 - start_layer->N->SE->W->end_layer, 1 - start_layer->E->SW->N->end_layer, ...
	//  4 ...  7 from end_layer to  start_layer: 4 - end_layer->N->SE->W->start_layer, 7 - end_layer->E->SW->N->start_layer, ...
	//  8 ... 15 - ortho-ortho-diagonal (turn right)
	//  8 ... 11 :  8 - start_layer->N->E->SW->end_layer,  9 - start_layer->E->S->NW->end_layer, ...
	// 12 ... 15 : 12 - end_layer->N->E->SW->start_layer, 13 - end_layer->E->S->NW->start_layer, ...
	// 16 ... 23 - ortho-ortho-diagonal (turn left)
	// 16 ... 19 : 16 - start_layer->N->W->SE->end_layer, 17 - start_layer->E->N->SW->end_layer, ...
	// 20 ... 23 : 20 - end_layer->N->W->SE->start_layer, 21 - end_layer->E->N->SW->start_layer, ...
	private void dbgCheck(){
//		if (nsTile== 630){
//			System.out.println("Conflict().dbgCheck, nsTile = "+nsTile);
//			System.out.println("Conflict().dbgCheck, nsTile = "+nsTile);
//		}
	}
	public Conflict(
			int nsTile,
			int start_layer,
			int end_layer,
			int start_dir) // for ortho-diag-ortho
	{
		this.nsTile = nsTile;
		if (end_layer > start_layer) {
			this.start_layer =               start_layer;
			this.end_layer =                 end_layer;
			this.start_dirs[start_dir] = true;
		} else {
			this.start_layer =               end_layer;
			this.end_layer =                 start_layer;
			this.start_dirs[start_dir + 4] = true;

		}
		dbgCheck();
	}

	Conflict(
			int nsTile,
			int start_layer,
			int end_layer,
			int start_dir,
			boolean right) // for ortho-ortho-diag
	{
		this.nsTile = nsTile;
		if (end_layer > start_layer) {
			this.start_layer =                 start_layer;
			this.end_layer =                   end_layer;
			this.start_dirs[(right?  8 : 16) + start_dir] = true;
		} else {
			this.start_layer =                 end_layer;
			this.end_layer =                   start_layer;
			this.start_dirs[(right? 12 : 20) + start_dir] = true;
		}
		dbgCheck();
	}

	public Conflict(int nsTile,
			int [] arr)
	{
		this.nsTile = nsTile;
		this.start_layer = arr[0];
		this.end_layer =   arr[1];
		for (int i = 0; i < start_dirs.length; i++){
			start_dirs[i] = (arr[2] & (1 << i)) != 0;
		}
	}
	Conflict(int bits)
	{
		this.start_layer = -1;
		this.end_layer =   -1;
		for (int i = 0; i < start_dirs.length; i++){
			start_dirs[i] = (bits & (1 << i)) != 0;
		}
		dbgCheck();
	}

	public int getSTile()
	{
		return this.nsTile;
	}

	@Override
	public String toString()
	{
		String s="Conflict"+
				" nsTile = "+getSTile()+
				" nl1 = "+ getStartLayer()+
				" nl2 = "+ getEndLayer()+
				" all = "+ getNumConflicts() + " ("+String.format("%06x", getDirBits())+")" +
				" odo = "+ getNumOrthoDiagOrthoConflicts() +
				" ood = "+ getNumOrthoOrthoDiagConflicts() +
				" number of odo incompatible triangles = "+ getIncompatibleOrthoDiagOrthoConflicts() +
				" number of odo dual triangles = " + getDualTriOrthoDiagOrthoConflicts();
		return s;
	}
	public boolean conflictExists(
			int start_dir4,
			boolean right,
			boolean reverse)
	{
		return start_dirs[8 + start_dir4 + (reverse? 4:0) + (right? 0 : 8)];
	}


	boolean combine (Conflict other_conflict)
	{
		dbgCheck();
		if ((other_conflict.start_layer == this.start_layer) && (other_conflict.end_layer == this.end_layer)) {
			for (int i = 0; i < start_dirs.length; i++) start_dirs[i] |= other_conflict.start_dirs[i];
			return true;
		}
		return false;
	}

	int getStartLayer(){
		return start_layer;
	}

	int getEndLayer(){
		return end_layer;
	}

	public int getStartLayer(boolean reverse){
		return reverse? end_layer: start_layer;
	}

	public int getEndLayer(boolean reverse){
		return reverse? start_layer : end_layer;
	}



	int getDirBits(){
		int dirs_bits = 0;
		for (int i = 0; i < start_dirs.length; i++){
			if (start_dirs[i]) dirs_bits |= (1 << i);
		}
		return dirs_bits;
	}

	int getDirBitsOrthoDiagOrtho(){
		int dirs_bits = 0;
		for (int i = 0; i < 8; i++){
			if (start_dirs[i]) dirs_bits |= (1 << i);
		}
		return dirs_bits;
	}

	int getDirBitsOrthoOrthoDiag(boolean right){
		int dirs_bits = 0;
		if (right) {
			for (int i = 0; i < 8; i++){
				if (start_dirs[i]) dirs_bits |= (1 << (i + 8));
			}
		} else {
			for (int i = 0; i < 8; i++){
				if (start_dirs[i]) dirs_bits |= (1 << (i + 16));
			}

		}
		return dirs_bits;
	}

	int [] toArray()
	{
		int dirs_bits = 0;
		for (int i = 0; i < start_dirs.length; i++){
			if (start_dirs[i]) dirs_bits |= (1 << i);
		}
		int [] rslt = {start_layer, end_layer, dirs_bits};
		return rslt;
	}

	int getNumConflicts(){
		int numbits = 0;
		for (int i = 0; i < start_dirs.length; i++){
			if (start_dirs[i]) numbits++;
		}
		return numbits;
	}

	int getNumOrthoDiagOrthoConflicts(){
		int numbits = 0;
		for (int i = 0; i < 8; i++){
			if (start_dirs[i]) numbits++;
		}
		return numbits;
	}

	int getNumOrthoOrthoDiagConflicts(){
		int numbits = 0;
		for (int i = 8; i < 24; i++){
			if (start_dirs[i]) numbits++;
		}
		return numbits;
	}

	int getIncompatibleOrthoDiagOrthoConflicts(){
		int num_incompat = 0;
		for (int i = 0; i < 8; i++){
			if (start_dirs[i]) {
				int [] incomp_bits= {
						i ^ 4,
						(i & 4) | ((i + 1) & 3),
						(i & 4) | ((i - 1) & 3)};
				for (int j = 1; j < 3; j++){ // skip dual triangles
					int i1 = incomp_bits[j];
					if (start_dirs[i1]){
						num_incompat ++;
					}
				}
			}
		}
		return num_incompat / 2;
	}

	public int getDualTriOrthoDiagOrthoConflicts()
	{
		int num_dual = 0;
		for (int i = 0; i < 4; i++) if (start_dirs[i] && start_dirs[i + 4]) num_dual++;
		return num_dual;
	}

	/**
	 * Get directions from center and layers of the 2 non-center tiles involved in each triangular conflict
	 * Full triangle is center:nl1 -> dir1 -> dir2 -> center:nl2
	 * @return array of 2-tuples:{dir1, dir2}
	 */
	public int [][] getInvolvedTiles(){
		int [][] involved = new int [getNumConflicts()][2];
		int nt = 0;
		for (int i = 0; i < 4; i++){
			if (start_dirs[i + 0]) {
				involved[nt][0] =   i * 2;
				involved[nt++][1] = ((i * 2) + 2) % 8;
			}
			if (start_dirs[i + 4]) {
				involved[nt][0] =   ((i * 2) + 2) % 8;
				involved[nt++][1] = i * 2;
			}
			if (start_dirs[i + 8]) {
				involved[nt][0] =   i * 2;
				involved[nt++][1] = ((i * 2) + 1) % 8;
			}
			if (start_dirs[i + 12]) {
				involved[nt][0] =   ((i * 2) + 1) % 8;
				involved[nt++][1] = i * 2;
			}

			if (start_dirs[i + 16]) {
				involved[nt][0] =   i * 2;
				involved[nt++][1] = ((i * 2) + 7) % 8;
			}
			if (start_dirs[i + 20]) {
				involved[nt][0] =   ((i * 2) + 7) % 8;
				involved[nt++][1] = i * 2;
			}
		}
		return involved;
	}

}
