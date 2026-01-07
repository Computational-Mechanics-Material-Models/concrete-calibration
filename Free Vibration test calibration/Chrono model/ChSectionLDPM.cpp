// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Erol Lale
// =============================================================================
// Section class for LDPM and CSL elements 
//
//  i)   Material
//  ii)  Projected area of facets
//  iii) Center of facets  
//  iv)  System of reference of facets
//
// A description of the material parameter can be found in: https://doi.org/10.1016/j.cemconcomp.2011.02.011
// =============================================================================

#include "ChSectionLDPM.h"

namespace chrono {
namespace fea {

// Construct CSL or LDPM section.

ChSectionLDPM::ChSectionLDPM( std::shared_ptr<ChMaterialVECT> material,  // material 
                                    double area,    // Projected total area of the facet
                                    ChVector<double> center,    // Center point of the facet area      
                                    ChMatrix33<double> facetFrame    /// local system of frame of facet 
                                       )
    : m_material{material}, m_area{area}, m_center{center}, m_facetFrame{facetFrame} {
        
}

ChSectionLDPM::ChSectionLDPM() {};

ChSectionLDPM::~ChSectionLDPM() {};


}  // end of namespace fea
}  // end of namespace chrono
