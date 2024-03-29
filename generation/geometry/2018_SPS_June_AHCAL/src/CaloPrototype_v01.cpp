//====================================================================
//  DD4hep Geometry driver for Sampling Calo BOX prototype
//--------------------------------------------------------------------
//  S.Lu, DESY
//  $Id:  $
//====================================================================
#include "DD4hep/Printout.h"
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"
#include "DDSegmentation/TiledLayerGridXY.h"
#include "LcgeoExceptions.h"

#include <iostream>
#include <vector>

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;
using namespace lcgeo ;

// workaround for DD4hep v00-14 (and older) 
#ifndef DD4HEP_VERSION_GE
#define DD4HEP_VERSION_GE(a,b) 0 
#endif


static Ref_t create_detector(LCDD& lcdd, xml_h element, SensitiveDetector sens)  {

  xml_det_t   x_det       = element;
  string      det_name    = x_det.nameStr();
  DetElement  sdet( det_name,x_det.id() );

  Layering    layering(x_det);

  // --- create an envelope volume and position it into the world ---------------------
  
  Volume envelope = XML::createPlacedEnvelope( lcdd,  element , sdet ) ;
  
  XML::setDetectorTypeFlag( element, sdet ) ;
  
  if( lcdd.buildType() == BUILD_ENVELOPE ) return sdet ;

  //-----------------------------------------------------------------------------------

  Material      air               = lcdd.air();

  sens.setType("calorimeter");

//====================================================================
//
// Read all the constant from compact.xml, user can update the value.
// Use them to build a calo box prototye.
//
//====================================================================

  double      Calo_half_x          = lcdd.constant<double>("HCAL_dim_x")/2.0;
  double      Calo_half_y          = lcdd.constant<double>("HCAL_dim_y")/2.0;
  double      Calo_half_z          = lcdd.constant<double>("HCAL_dim_z")/2.0;
  double      Calo_Layer_ncell_x   = lcdd.constant<int>("HCAL_nCells_x");
  double      Calo_Layer_ncell_y   = lcdd.constant<int>("HCAL_nCells_y");

  printout( DD4hep::DEBUG,  "building SamplingCaloBoxPrototype_v01",
	    "Calo_half_x : %e    Calo_half_y: %e    Calo_half_z: %e ",
  	    Calo_half_x, Calo_half_y, Calo_half_z) ;


//====================================================================
//
// general calculated parameters
//
//====================================================================
 

//====================================================================
//
// Chambers in the CaloBox
//
//====================================================================

    int layer_num = 0;
    int layerType = 0;
    double cal_hx = 0;
    double cal_hy = 0;
    double cal_hz = 0;

    double layer_pos_z = - Calo_half_z;

    for (xml_coll_t c(x_det, _U(layer)); c; ++c) {
        xml_comp_t x_layer = c;
        int repeat = x_layer.repeat();                // Get number of times to repeat this layer.
        const Layer* lay = layering.layer(layer_num); // Get the layer from the layering engine.
        double layer_thickness = lay->thickness();

	//int layerType  = lay->layerType();
        string layer_type_name   = _toString(layerType,"layerType%d");
	
	//calorimeter dimensions
	cal_hx = (double) (Calo_Layer_ncell_x * 3.0)/2.;
	cal_hy = (double) (Calo_Layer_ncell_y * 3.0)/2.;

        // Loop over repeats for this layer.
        for (int j = 0; j < repeat; j++) {
            string layer_name = _toString(layer_num, "layer%d");
            DetElement layer(layer_name, layer_num);

            Volume layer_vol(layer_type_name, Box(cal_hx, cal_hy, layer_thickness / 2), air);
            
            // Create the slices (sublayers) within the layer.
            double slice_pos_z = -(layer_thickness / 2);
            int slice_number = 0;
            
            for (xml_coll_t k(x_layer, _U(slice)); k; ++k) {
                xml_comp_t x_slice = k;
                string slice_name = _toString(slice_number, "slice%d");
                double slice_thickness = x_slice.thickness();
                Material slice_material = lcdd.material(x_slice.materialStr());

                slice_pos_z += slice_thickness / 2;
                // Slice volume & box
                Volume slice_vol(slice_name, Box(cal_hx, cal_hy, slice_thickness / 2), slice_material);
                
                
                if (x_slice.isSensitive()) {
                    sens.setType("calorimeter");
                    slice_vol.setSensitiveDetector(sens);                   
                } 

                // Set region, limitset, and vis.
                slice_vol.setAttributes(lcdd, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
                // slice PlacedVolume
                layer_vol.placeVolume(slice_vol, Position(0, 0, slice_pos_z));
                
                // Increment Z position for next slice.
                slice_pos_z += slice_thickness / 2;
                // Increment slice number.
                ++slice_number;
            }
            
            
            // Set region, limitset, and vis.
            layer_vol.setAttributes(lcdd, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());
           
            // Layer position in Z within the stave.
            layer_pos_z += layer_thickness / 2;
            // Layer physical volume.
            PlacedVolume layer_phv = envelope.placeVolume(layer_vol, Position(0, 0, layer_pos_z));
            //layer_phv.addPhysVolID("layer", layer_num);
            layer_phv.addPhysVolID("K", layer_num);
            layer.setPlacement(layer_phv);
            
            // Increment the layer Z position.
            layer_pos_z += layer_thickness / 2;
            // Increment the layer number.
            ++layer_num;
	++layerType;

        }

    }

 
    return sdet;

}

DECLARE_DETELEMENT(CaloPrototype_v01, create_detector)
