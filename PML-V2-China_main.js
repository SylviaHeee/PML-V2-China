/**
 * How to cite. 
 * @reference
 * 1. He, S., Zhang, Y., Ma, N., Tian, J., Kong, D., and Liu, C.: 
 *    A daily and 500 m coupled evapotranspiration and gross primary production product across China during 2000–2020, 
 *    Earth Syst. Sci. Data Discuss. [preprint], https://doi.org/10.5194/essd-2022-183, in review, 2022.
 * 
**/

var pkg_PML = {};
var dataset = {};  

var I_interp = true;
var meth_interp = 'bilinear'; // or 'bicubic'; for meteometeorological forcing spatial interpolatation

function init_dataset() {
    if (!pkg_main.is_empty_dict(dataset)) return;
    // var I_interp = true; // whether Interpolate MODIS LAI, Emissivity and Albedo
    // `meth_interp` is used to resample  into high-resolution
    // not suggest 'biculic'. bicubic can't constrain values in reasonable boundary.
    
    var filter_date_all = ee.Filter.date('2000-02-26', '2021-01-01');

   /** fix MCD12Q1_006 land cover code. */
    var ImgCol_land = imgcol_land.select(0).map(function (land) {
        //for MCD12Q1_006 water and unc type is inverse
        land = land.remap([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
            [17, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 0]);
        return (land);
    });

    /** 1.1 CMFD and CO2 */
    var ImgCol_co2 = co2.toList(co2.size())
        .map(function (f) {
            f = ee.Feature(f);
            var date = ee.Date.parse('YYYY-MM-dd', f.get('date'));
            // print(date);
            return ee.Image.constant(f.get('average'))
                .toFloat()
                .set('system:time_start', date.millis())
                .set('system:id', date.format('YYYY-MM-dd'))
                .set('system:index', date.format('YYYY-MM-dd'));
        });
    ImgCol_co2 = ee.ImageCollection(ImgCol_co2).select([0], ['co2'])
        .filter(filter_date_all)
        .sort("system:time_start");
    // print(ImgCol_co2)
      
    // **China Or Globe;i.e., using CMFD or GLDAS**
    Imgcol_cmfd = Imgcol_cmfd.merge(GLDAS_daily)
                  .select(['Tavg','q','U2','Prcp','pa','Rs','Rln','Tmax','Tmin'])
                  .filter(filter_date_all);
    
    // **ERA5-LST**
    /**bands2imgcol for LST
     * multiple bands image convert to image list
     * 
     * The bandName should be like that "b20020101_LST".
     * 
     * @param  {[type]} img      multiple bands image
     * @param  {[type]} bandname the new bandname
     * @return {ee.List}         List of images
     */
    var bands2imgcol_lst = function(img, bandname) {
        bandname = bandname || "b1";
        img = ee.Image(img);
        var names = img.bandNames(); //ee.List()
        var n = names.size();
    
        var imgcol_lst = names.map(function (name) {
            var date = ee.Date.parse('YYYYMMdd', ee.String(name).slice(1, 9));
            return img.select([name], [bandname])
                .set('system:time_start', date.millis())
                .set('system:id', date.format('yyyy_MM_dd'))
                .set('system:index', date.format('yyyy_MM_dd'));
        });
        return imgcol_lst;
    }
    var ImgCol_LST = ee.ImageCollection(Daily_LST.toList(30).map(function(img){ 
        return bands2imgcol_lst(img, 'LST');
    }).flatten()).filter(filter_date_all)
    
    /** 1.2 MODIS products: LAI, Albedo, Emissivity  */
    var imgcol_lai;
    if (I_interp) {
       imgcol_lai = ee.ImageCollection(newLAI.toList(2).map(function(img){
                          return pkg_main.bands2imgcol(img, 'LAI');
                        }).flatten());

        imgcol_emiss = ee.ImageCollection(imgcol_emiss.toList(1000))
            .map(function (img) {
                var emiss = img.select(0).expression('b() * 0.002 + 0.49'); //.toFloat(); //.toUint8()
                return img.select('qc').addBands(emiss);
            }).select([1, 0], ['Emiss', 'qc']);

        imgcol_albedo = ee.ImageCollection(imgcol_albedo.toList(1000))
            .map(function (img) {
                var albedo = img.select(0).multiply(0.001); //.toFloat();
                return img.select(1).addBands(albedo);
            }).select([1, 0], ['Albedo', 'qc']);//scale factor 0.001, no units;
        // print('Interped');
    } else {
        /** No Interpolation MODIS INPUTS */
        imgcol_lai = ee.ImageCollection('MODIS/006/MCD15A3H').select('Lai')
            .map(function (img) { return img.multiply(0.1).copyProperties(img, img.propertyNames()); }); //scale factor 0.1

        imgcol_emiss = ee.ImageCollection('MODIS/006/MOD11A2')
            .select(['Emis_31', 'Emis_32'])
            .map(function (img) {
                return img.reduce(ee.Reducer.mean()).multiply(0.002).add(0.49)
                    .copyProperties(img, ['system:time_start', 'system:id']);
            }).select([0], ['Emiss']);

        var Albedo_raw = ee.ImageCollection('MODIS/006/MCD43A3').select(['Albedo_WSA_shortwave'])
            .map(pkg_trend.add_dn(true));
        imgcol_albedo = pkg_agg.aggregate_prop(Albedo_raw, 'd8', 'mean')
            .map(function (img) { return img.addBands(img.multiply(0.001)).select([1]); })
            .select([0], ['Albedo']);
        // print('No Interped');
    }
    
    dataset = {
        ImgCol_land  : ImgCol_land, 
        Imgcol_cmfd  : Imgcol_cmfd,
        ImgCol_LST   : ImgCol_LST,
        ImgCol_co2   : ImgCol_co2, 
        imgcol_lai   : imgcol_lai,
        imgcol_emiss : imgcol_emiss, 
        imgcol_albedo: imgcol_albedo,
    };
    // print(dataset)
}

/** 1. fill_options --------------------------------------------------------- */
// default_true(null)
// default_true(undefined)
function default_true(x) {
    return (x === undefined || x === null) ? true : x;
}

var prj = pkg_export.getProj(PML_V2);
var options = {
    Export: {
        range: [73,18,136,54],//yearly [69,15,140,55],//[-180, -60, 180, 89]world,
        cellsize: 500,// meaningless   1 / 240, //1/240,
        type: 'asset',  //drive    //asset
        crs: 'EPSG:4326',//'SR-ORG:6974', //projects/pml_evapotranspiration
        crsTransform: prj.crsTransform    
    }
};
// print(options);

/**
 * [fill_options description]
 *
 * @param {Dictionary} opt
 * - `timescale`: `daily` or `yearly`
 * - `year_begin`:
 * - `year_end`:
 * - `is_PMLV2`:
 * - `is_save`:
 * - `is_dynamic_lc`:
 * - `timescale`:
 * - `folder`:
 * @return {[type]} [description]
 *
 * @example
 var opt = {
    timescale = 'yearly',
    folder = 'projects/pml_evapotranspiration/PML/V2/yearly'
 }
 */
pkg_PML.fill_options = function(opt, verbose){

    var year_begin    = opt.year_begin || 2000;//2000;
    var year_end      = opt.year_end   || 2020;
    var is_PMLV2      = default_true(opt.is_PMLV2);
    var is_save       = default_true(opt.is_save);
    var is_dynamic_lc = default_true(opt.is_dynamic_lc);
    var timescale     = opt.timescale || "yearly";

    var bands, folder2;

    if (is_PMLV2) {
        bands = ['GPP', 'Ec', 'Es', 'Ei', 'ET_water', 'qc']; //,'qc'
    } else {
        bands = ['Ec', 'Es', 'Ei', 'ET_water', 'qc'];
    }
    folder2 = opt.folder; // if folder not null, than replace
    options.Export.folder = folder2;

    var prefix = is_PMLV2 ? "PMLV2_" : "PMLV1_";
    prefix = prefix.concat(timescale).concat("_");
    options.prefix = prefix;
  
    options.is_save       = is_save;
    options.year_begin    = year_begin; 
    options.year_end      = year_end; 
    options.is_PMLV2      = is_PMLV2; 
    options.is_dynamic_lc = is_dynamic_lc;    
    options.timescale     = timescale; 
    options.bands         = bands;   
    
    if (verbose) print('PML options', options);
};

/** 2. ------------------------------------------------------------------- */
/**
 * Prepare INPUT datset for PML_V2
 *
 * @param {[type]} begin_year [description]
 * @param {[type]} end_year   [description]
 */

/** PML GLOBAL PARAMETERS */
var Gsc = 0.0820,  // solar constant in unit MJ m-2 min-1,
    as = 0.25,    // parameter Rs/Ra=as+bs*n/N; calibration from our solar radiation measurement
    bs = 0.50,    // parameter Rs/Ra=as+bs*n/N;
    alfa = 0.23,    // surface albedo of grass
    alfa_forest = 0.22,    // surface albedo of forest
    alfa_crop = 0.14,    // surface albedo of crop

    kmar = 0.40,    // von Karman's constant 0.40 
    Zob = 15,      // m, making sure higher than hc
    Cp = 1.0164,  // 4.2 * 0.242, specific heat at constant pressure, 1.013  [J g-1 0C-1]
    epsl = 0.622,   // ratio molecular weight of water vapour/dry air

    /** PML_v1 parameters for Gc */
    kQ = 0.4488,  // extinction coefficient
    kA = 0.7,     // the attenuation of net all-wave irradicance, typically about 0.6-0.8 (Denmend, 1976, Kelliher FM et al., (1995))
    Q50 = 30,      // the value of absorbed PAR when gs=gsx/2, W/m2
    D0 = 0.7;     // the value of VPD when stomtal conductance is reduced  kpa 

/**
 * SEVEN OPTIMIZED PARAMETERS
 * 
 * Alpha  : initial photochemical efficiency, 0.02-0.08
 * Thelta : the initla slope of the slope of CO2 response curve[umol m-2 s-1]/[umol mol-1], 1
 * m      : Ball-Berry coefficient 2-20
 * Am_25  : the maximum catalytic capacity of Rubisco per unit leaf area at 25 degree
 * kQ     : the value of VPD when stomtal conductance is reduced 
 * kA     : extinction coefficient
 *
 * TWO INTERCEPTION PARAMETERS
 * S_sls  : specific canopy rainfall storage capacity per unit leaf area (mm)
 * fER0   : 
set:
13 (Urban and Built-Up)           = 5  (mixed forest)
16 (Barren or Sparsely Vegetated) = 10 (grassland)
 */
 
/** canopy height */

var hc_raw = ee.List([0.01, 10, 15, 10, 15, 10,
    1.5, 1.5, 5, 5, 0.2, 1,
    0.5, 10, 1, 0.01, 0.05, 0.1]); //update  2022/3/18


//update  2022/3/28
var Alpha_raw = ee.List([
0,	0.0372233268591397,	0.0388628012500654,	0.0372233268591397,	0.0388020938048238,	0.0388020938048238,
0.0289228534310204,	0.0289228534310204,	0.0388020938048238,	0.0392246392752068,	0.0499410977946029,
0.0280720284215188,	0.0391362599023089,	0.029,	0.0391362599023089,	0,	0.012494558253439,	0
]);

var Thelta_raw = ee.List([
0,	0.0429386162358163,	0.00924519672646321,	0.0429386162358163,	0.0223237210297614,	0.0223237210297614,
0.0627125368170904,	0.0627125368170904,	0.0223237210297614,	0.0150555543993891,	0.0686687007401534,
0.0311597420693756,	0.0599394333108481,	0.036,	0.0599394333108481,	0,	0.0415726727694529,	0
]);

var m_raw = ee.List([
0,	3.7258703861433,	7.99962937514182,	3.7258703861433,	7.79915432528032,	7.79915432528032,
6.31798010167713,	6.31798010167713,	7.79915432528032,	6.36436208511253,	13.7066752927126,
21.8661000019827,	4.97184217133699,	8.539,	4.97184217133699,	0,	5.53587362716509,	0
]);

var Am_25_raw = ee.List([
0,	46.3662066113467,	9.65501396148521,	46.3662066113467,	9.65515714791474,	9.65515714791474,
12.8114224918777,	12.8114224918777,	9.65515714791474,	2.39920334368553,	6.82002235382988,
46.3697409005725,	29.9988953534116,	15.93,	29.9988953534116,	0,	46.3572383031857,	0
]);


var D0_raw = ee.List([
0.7,	1.98700140635691,	0.977880421746307,	1.98700140635691,	0.79920463685339,	0.79920463685339,
0.898710865055111,	0.898710865055111,	0.79920463685339,	0.882831813769275,	0.51611202868386,
1.60425965286808,	1.9924289725051,	0.501,	1.9924289725051,	0.7,	1.9952052290022,	0.7
]);


var kQ_raw = ee.List([
0.6,	0.999991726748929,	0.54474515567827,	0.999991726748929,	0.505804693924446,	0.505804693924446,
0.536203249602075,	0.536203249602075,	0.505804693924446,	0.423092694367454,	0.998064378036293,
0.105410133969185,	0.203023793715866,	0.593,	0.203023793715866,	0.6,	0.995399825050184,	0.6
]);


var kA_raw = ee.List([
0.7,	0.682864390304589,	0.800212962972644,	0.682864390304589,	0.887793586008346,	0.887793586008346,
0.142592616336453,	0.142592616336453,	0.887793586008346,	0.889024191280602,	0.897151751732363,
0.889996569931878,	0.885645491530372,	0.68,	0.885645491530372,	0.7,	0.693649385520565,	0.7
 
]);


var S_sls_raw = ee.List([
0,	0.167414335926914,	0.115517312854426,	0.167414335926914,	0.0586943663801799,	0.0586943663801799,
0.159778248924219,	0.159778248924219,	0.0586943663801799,	0.0546372902881614,	0.127605412582778,
0.000511268729862606,	0.00906152882973127,	0.131,	0.00906152882973127,	0,	0.169720305292886,	0
]);


var fER0_raw = ee.List([
0,	0.00735002832672545,	0.0156231003452271,	0.00735002832672545,	0.0171455985347659,	0.0171455985347659,
0.0696392060761715,	0.0696392060761715,	0.0171455985347659,	0.145912246955316,	0.00306375320798436,
0.00537734487657938,	0.0105707049982767,	0.01,	0.0105707049982767,	0,	0.146856871618889,	0
]);


var VPDmin_raw = ee.List([
1,	0.650279317925201,	0.650072934315015,	0.650279317925201,	0.730741647829743,	0.730741647829743,
0.789127079781859,	0.789127079781859,	0.730741647829743,	1.41032689458432,	1.48462591306493,
0.651729624915632,	1.39361527264971,	0.657,	1.39361527264971,	1,	1.48946604240688,	1
]);


var VPDmax_raw = ee.List([
4,	5.32480342745534,	4.93453439271283,	5.32480342745534,	5.4764336258539,	5.4764336258539,
3.50212148818488,	3.50212148818488,	5.4764336258539,	6.4997584702857,	6.47208792098744,
5.51369884363017,	4.99729923184796,	3.5,	4.99729923184796,	4,	6.4962560535831,	4
]);
/**
 * Construct parameters depend on landcover type
 *
 * @param  {ee.Image} landcover [description]
 * @param  {ee.List}  list      [description]
 * @return {ee.Image}           [description]
 */
function propertyByLand_v2(landcover, list) {
    landcover = ee.Image(landcover);
    // modis landcover 18 types
    var lands = ee.List.sequence(0, 17).map(function (i) {
        i = ee.Number(i);
        var land = landcover.eq(i).float();
        var prop = ee.Number(list.get(i));
        return land.multiply(prop);
    });
    return ee.ImageCollection(lands).sum();
}

/** Vapor Pressure in kPa with temperature in degC */
function vapor_pressure(t) {
    return t.expression('0.6108 * exp(17.27 * b() / (b() + 237.3))');
}

/**
 * PML_V2 (Penman-Monteith-Leuning) model
 *
 * sub functions:
 *     -- PML_daily(img)
 *     `PML_daily` has all the access of yearly land cover based parameters 
 *     (e.g. gsx, hc, LAIref, S_sls). 
 *     
 * 
 *     -- PML_year(INPUTS)
 *     
 * @param {Integer} year Used to filter landcover data and set landcover depend parameters.
 * @param {boolean} is_PMLV2 Default is true, and PML_V2 will be used. If false, 
 *                     PML_V1 will be used.
 * 
 * @return {ee.ImageCollection} An ImageCollection with the bands of 
 *                                 ['GPP', 'Ec', 'Es', 'Ei', 'ET_water','qc'] for PML_V2;
 *                                 ['Ec', 'Es', 'Ei', 'ET_water','qc'] for PML_V1;
 * 
 */
function PML(year, is_PMLV2) {
    var year_max = 2020,
        year_min = 2001;
    var year_land = year;
    if (year > year_max) year_land = year_max;
    if (year < year_min) year_land = year_min;

    var filter_date_land = ee.Filter.calendarRange(year_land, year_land, 'year');
    var land = ee.Image(dataset.ImgCol_land.filter(filter_date_land).first()).aside(print,"land")//.aside(Map.addLayer,{},'China_500m',false); //land_raw was MODIS/051/MCD12Q1

    /** remove water, snow and ice, and unclassified land cover using updateMask */
    // var mask     = land.expression('b() != 0 && b() != 15 && b() != 17');
    // land         = land.updateMask(mask);
    // var landmask = ee.Image(1).updateMask(mask);

    /** Initial parameters */
    // gsx, hc, LAIref, S_sls can be accessed by `PML_daily`, in its parent env
    var hc = propertyByLand_v2(land, hc_raw);

    if (is_PMLV2) {
        var Alpha = propertyByLand_v2(land, Alpha_raw),
            Thelta = propertyByLand_v2(land, Thelta_raw),
            m  = propertyByLand_v2(land, m_raw),
            Am_25 = propertyByLand_v2(land, Am_25_raw);
        // Ca      = 380; //umol mol-1
        D0 = propertyByLand_v2(land, D0_raw);
        kQ = propertyByLand_v2(land, kQ_raw);
        kA = propertyByLand_v2(land, kA_raw);
        var VPDmin = propertyByLand_v2(land, VPDmin_raw);
        var VPDmax = propertyByLand_v2(land, VPDmax_raw);
        // VPDmin = 0.93; VPDmax = 4.3;
        // for PML_v1 D0, kQ, kA are constant parameters.
    }else{
      var gsx = propertyByLand_v2(land, gsx_raw)   //only for PML_v1
    }
    
    // parameters for Ei
    var LAIref = ee.Image(5), //propertyByLand_v2(land, LAIref_raw),
        S_sls = propertyByLand_v2(land, S_sls_raw),
        fER0 = propertyByLand_v2(land, fER0_raw);

    /**
    * Calculate daily PML GPP and ET using CMFD and MODIS inputs.
    * 
    * @param  {Image} img CMFD meteorological forcing data and MODIS remote sensing data
    *    with bands: ['LAI', 'Emiss', 'Albedo', 'pa', 'Tmax', 'Tmin', 'Tavg', 'Prcp', 'Rln', 'Rs', 'U2']
    * 
    * @return {Image} PML_ET with bands of ['ET_water', 'Es_eq', 'Ec', 'Ei', 'Pi']; 
    *                 If v2 = true, GPP also will be returned.
    */
    function PML_daily(img) {
        img = ee.Image(img);
        var Ca = img.select('co2');   //umol mol-1
        var q = img.select('q');     // kg/kg;
        var p = img.select('pa');    // kPa
        var u2 = img.select('U2');    // m/s

        var Tmax = img.select('Tmax');  // degC
        var Tmin = img.select('Tmin');  // degC
        var Tavg = img.select('Tavg');  // degC
        var LST = img.select('LST');   // degC

        var Rln = img.select('Rln');   // W/m2/s, not MJ/m2/d 
        var Rs = img.select('Rs');    // W/m2/s

        var albedo = img.select('Albedo');// %
        var emiss = img.select('Emiss'); // %
        var LAI = img.select('LAI');   // 0 - 

        // var lambda = 2500; // latent heat of vaporization, 2500 [J g-1]  at 25 degC //KDD
        // lambda = Tavg.multiply(-2.2).add(lambda);//KDD
        var lambda = 2501;//Shaoyang
        lambda = Tavg.multiply(-2.361).add(lambda);//Shaoyang
        
        /** 
        * ACTUAL VAPOUR PRESSURE
        * https://www.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html, Eq-17
        */
        var ea = img.expression('q * p / (0.622 + 0.378 * q)', { 'p': p, 'q': q });

        // saturation vapour pressure from Tair
        var es_tmax = vapor_pressure(Tmax);
        var es_tmin = vapor_pressure(Tmin);
        var es_tavg = vapor_pressure(Tavg);
        var es = es_tmax.add(es_tmin).divide(2);

        var VPD = es.subtract(ea).max(0.001);

        var rou_a = img.expression('3846 * Pa / (Tavg + 273.15)',
            { 'Pa': p, 'Tavg': Tavg });
        var gama = img.expression('Cp*Pa/(0.622*lambda)',
            { Cp: Cp, Pa: p, lambda: lambda }); // kpa/0C
        
        var slop = img.expression('4098 * es_tavg / pow(Tavg + 237.3, 2)',
            { 'es_tavg': es_tavg, 'Tavg': Tavg });
        
        // downward Solar Radiation
        var Stefan = 4.903e-9;// Stefan-Boltzmann constant [MJ K-4 m-2 day-1],
        var Rns = ee.Image(1).subtract(albedo).multiply(Rs);
        
        // if no LST, using Tavg 
        // var RLout = img.expression('Emiss * Stefan * pow(Tavg+273.15, 4)',
        //     { 'Emiss': emiss, Stefan: Stefan, Tavg: Tavg }).divide(0.0864);
        
        var RLout = img.expression('Emiss * Stefan * pow(LST+273.15, 4)',
            { 'Emiss': emiss, Stefan: Stefan, LST: LST }).divide(0.0864);
            
        // var Rnl = Rln.subtract(RLout);//KDD
        var Rnl = emiss.multiply(Rln).subtract(RLout);//Shaoyang
        var Rn = Rns.add(Rnl).max(0.0);    // to ensure Rn >= 0;
        var PAR = Rs.multiply(0.45).max(0); // could be used modis data to replace
        // units convert: http://www.egc.com/useful_info_lighting.php

        var Gc, GPP;
        var fvpd_gc = VPD.expression('1/(1+b()/D0)', { D0: D0 });        // leuning
        // var fvpd = VPD.expression('exp(-D0 * pow(b(), 2))', {D0:D0}); // yongqiang, f_VPD = exp(-D0 * VPD.^2);
        // var VPD_sqrt = VPD.sqrt();
        // var fvpd = VPD_sqrt.expression('b()*(b() < 1) + 1/b()*(b() >= 1)');
        var fvpd = VPD.expression('(VPDmax - b())/(VPDmax - VPDmin)', { VPDmin: VPDmin, VPDmax: VPDmax })
            .min(1.0).max(0.0);

        if (is_PMLV2) {
            var PAR_mol = PAR.multiply(4.57);    // from [W m-2] to [umol m-2 s-1]

            /** G flux part */
            var fT2 = Tavg.expression('exp(0.031*(b()-25))/(1 +exp(0.115*(b()-41)))').min(1.0);
            var Am = Am_25.multiply(fT2)//Shaoyang

            var P1 = Am.multiply(Alpha).multiply(Thelta).multiply(PAR_mol),
                P2 = Am.multiply(Alpha).multiply(PAR_mol),
                P3 = Am.multiply(Thelta).multiply(Ca),
                // P4 = Alpha.multiply(Thelta).multiply(PAR_mol).multiply(Ca).divide(fT2);//KDD
                P4 = Alpha.multiply(Thelta).multiply(PAR_mol).multiply(Ca);//Shaoyang

            var Ags = P1.expression('Ca*P1/(P2*kQ + P4*kQ) * (kQ*LAI + log((P2+P3+P4)/(P2+P3*exp(kQ*LAI) + P4)))', 
                { Ca: Ca, P1: P1, P2: P2, P3: P3, P4: P4, kQ: kQ, LAI: LAI, fT2: fT2 });  // umol cm-2 s-1
            GPP = Ags.multiply(1.0368).multiply(fvpd).rename('GPP'); //86400/1e6*12

            var img_check = GPP.addBands([rou_a, gama, slop, PAR, PAR_mol, fT2, P1, P2, P3, P4])
                .rename(['gpp', 'rou_a', 'gama', 'slop', 'par', 'par_mol', 'fT2', 'p1', 'p2', 'p3', 'p4']);

            // Gc = m.expression('m/Ca*Ags*1.6*fvpd_gc', { m: m, Ca: Ca, Ags: Ags, fvpd_gc: fvpd_gc });//!!!back Shaoyang
            Gc = m.expression('m/Ca*Ags*fvpd*1.6*fvpd_gc', { m: m, Ca: Ca, Ags: Ags, fvpd_gc: fvpd_gc,fvpd:fvpd });
            
            // Convert from mol m-2 s-1 to cm s-1 to m s-1
            Gc = Gc.expression('Gc*1e-2/(0.446*(273/(273+Tavg))*(Pa/101.3))',
                { Gc: Gc, Tavg: Tavg, Pa: p }); // unit convert to m s-1
        } 
        else {
            // Conductance and ET component
            Gc = LAI.expression('gsx/kQ*log((PAR+Q50)/(PAR*exp(-kQ*LAI)+Q50))*fvpd_gc',
                { gsx: gsx, kQ: kQ, PAR: PAR, Q50: Q50, LAI: LAI, fvpd_gc: fvpd_gc });
        }
        Gc = Gc.max(1e-6);
        // known bug: bare, ice & snow, unc, all zero parameters will lead to p1, p2, p3, p4 = 0,
        //            GPP = 0/0(masked), and Ec = masked.

        /** AERODYNAMIC CONDUCTANCE */
        // var d = hc.multiply(0.64);//KDD
        // var zom = hc.multiply(0.13);//KDD
        // var zoh = zom.multiply(0.1);//KDD
        var d = hc.multiply(0.667);//Ning
        var zom = hc.multiply(0.125);//Ning
        var zoh = zom.multiply(0.135);//Ning
        
        var uz = img.expression('log(67.8*Zob - 5.42)/4.87 * u2',
            { Zob: Zob, u2: u2 });
        var Ga = img.expression('uz*kmar*kmar / (log((Zob-d)/zom) * log((Zob-d)/zoh))',
            { uz: uz, kmar: kmar, Zob: Zob, zom: zom, zoh: zoh, d: d });

        // Equilibrium evaporation
        var Eeq = img.expression('slop/(slop+gama)*Rn', { slop: slop, gama: gama, Rn: Rn })
            .divide(lambda).multiply(86.4) // convert W/m2/s into mm
            .max(0.0001);
        // Penman Monteith potential ET
        var Evp = VPD.expression('(gama/(slop+gama))*((6430 * (1 + 0.536*u2) * VPD)/lambda)',
            { slop: slop, gama: gama, u2: u2, VPD: VPD, lambda: lambda })
            .max(0);
        var mask_water = land.expression('b() == 0 || b() == 15'); //water, snow&ice
        var ET_water = Eeq.add(Evp).updateMask(mask_water).rename('ET_water');

        // // Convert MJ/m2/day into W/m2;
        // Rn  = Rn.divide(0.0864).max(0);
        // PAR = PAR.divide(0.0864).max(0);

        // Conductance and ET component
        var Tou = LAI.expression('exp(-kA*LAI)', { kA: kA, LAI: LAI });

        // % Transpiration from plant cause by radiation water transfer
        var LEcr = slop.expression('slop/gama*Rn *(1 - Tou)/(slop/gama + 1 + Ga/Gc)',
            { slop: slop, gama: gama, Rn: Rn, Tou: Tou, Ga: Ga, Gc: Gc });               // W/m2
        // var LEcr = landmask.* LEcr;

        // % Transpiration from plant cause by aerodynamic water transfer
        var LEca = slop.expression('(rou_a * Cp * Ga * VPD / gama)/(slop/gama + 1 + Ga/Gc)',
            { rou_a: rou_a, Cp: Cp, Ga: Ga, Gc: Gc, VPD: VPD, gama: gama, slop: slop }); // W/m2

        // % making sure vegetation transpiration is negaligable, this is very important for very dry Sahara
        // Should take it seriously. LAI = 0, will lead to a extremely large value. 
        // Update 24 Aug'2017, kongdd
        LEca = LEca.where(LAI.lte(0.0), 0.0);
        LEcr = LEcr.where(LAI.lte(0.0), 0.0);
        var LEc = LEca.add(LEcr);

        // % Soil evaporation at equilibrium
        var LEs_eq = slop.expression('(slop/gama)* Rn *Tou/(slop/gama + 1)',
            { slop: slop, gama: gama, Rn: Rn, Tou: Tou });

        /** W/m2 change to mm d -1 */
        var coef_MJ2mm = lambda.divide(86.4); // ET./lambda*86400*10^-3;
        var Es_eq = LEs_eq.divide(coef_MJ2mm);
        var Ecr = LEcr.divide(coef_MJ2mm);
        var Eca = LEca.divide(coef_MJ2mm);
        var Ec = LEc.divide(coef_MJ2mm);

        /** 
        * Interception Precipitation Evaporation: prcp_real = prcp - Ei 
        * @references 
        * Van Dijk, A.I.J.M. and Warren, G., 2010. The Australian water resources assessment system. Version 0.5, 3(5). P39
        */
        var prcp = img.select('Prcp');
        var fveg = LAI.expression('1 - exp(-LAI/LAIref)', { LAI: LAI, LAIref: LAIref });
        var Sveg = S_sls.multiply(LAI);
 
        var fER = fveg.multiply(fER0);
        var prcp_wet = LAI.expression('-log(1 - fER0) / fER0 * Sveg / fveg',
            { fER0: fER0, fveg: fveg, Sveg: Sveg });
        var Ei = LAI.expression('(P < Pwet) * fveg * P + (P >= Pwet) * ( fveg*Pwet + fER*(P - Pwet) )',
            { fveg: fveg, fER: fER, P: prcp, Pwet: prcp_wet });
        var Pi = prcp.subtract(Ei);
        // (P < Pwet) * fveg * P + (P >= Pwet) * ( fveg*Pwet + fER*(P - Pwet) )
        //    NA and infinite values should be replaced as zero. But GEE where and 
        //    updatemask are incompetent.
        // ----------------------------------------------------------------------

        // var newBands = ['ETsim', 'Es', 'Eca', 'Ecr', 'Ei', 'Eeq', 'Evp', 'Es_eq'];
        var newBands = ['Es_eq', 'Ec', 'Ei', 'Pi']; //'Eeq', 'Evp', 'ETsim', 'Es'
        var newImg = ee.Image([Es_eq, Ec, Ei, Pi]).rename(newBands);
        if (is_PMLV2) newImg = newImg.addBands(GPP); //PML_V2

        newImg = newImg.updateMask(mask_water.not()).addBands(ET_water); //add ET_water
        // Comment 2018-09-05, to get yearly sum, it can be converted to uint16
        // otherwise, it will be out of range.
        newImg = newImg.multiply(1e2).toUint16(); //CONVERT INTO UINT16 

        if (I_interp) {
            var qc = img.select('qc');
            newImg = newImg.addBands(qc);
        }

        var beginDate = ee.Date(img.get('system:time_start'));
        return pkg_main.setImgProperties(newImg, beginDate);
        // return pkg_main.setImgProperties(img_check, beginDate);
    }

    /**
    * Calculate a period PML
    *
    * @param {ee.ImageCollection} INPUTS Multibands ImageCollection returned 
    * by PML_INPUTS_daily
    */
    
    function PML_INPUTS_daily(begin_year, end_year) {
      if (typeof end_year === 'undefined') { end_year = begin_year; }
      begin_year = ee.Number(begin_year);
      end_year = ee.Number(end_year).add(ee.Number(1));
  
      var begin_yearStr = begin_year.format('%d'),
          end_yearStr = end_year.format('%d');
          
      var date_begin = ee.Date(ee.Algorithms.If(begin_year.eq(ee.Number(2000)),
          begin_yearStr.cat("-02-26"), begin_yearStr.cat("-01-01"))),
          date_end = ee.Date(end_yearStr.cat("-01-01"));
      
      var filter_date = ee.Filter.date(date_begin, date_end);
      
      // **MODIS**
      var LAI_d8 = dataset.imgcol_lai.filter(filter_date);//.merge(lai_miss);
      var Albedo_d8 = dataset.imgcol_albedo.filter(filter_date);
      var Emiss_d8 = dataset.imgcol_emiss.filter(filter_date);
      var modis_input = pkg_join.SaveBest(Emiss_d8, LAI_d8);
      modis_input = pkg_join.SaveBest(modis_input, Albedo_d8);
      // print(modis_input);
      
      if (I_interp) {
          // add qc bands
          modis_input = modis_input.map(function (img) {
              var qc = img.expression('b("qc") + b("qc_1")*8').toUint8(); //qc, 0-2:emiss, 3-5:albedo
              return img.select(['LAI', 'Emiss', 'Albedo']).addBands(qc);
          });
      }
      
      //** CO2 **
      var CO2_d8 = ee.ImageCollection(dataset.ImgCol_co2)
                  .filter(filter_date)
                  .select([0], ['co2'])
                  .sort("system:index");
      
      var modis_input_co2 = pkg_join.SaveBest(modis_input, CO2_d8);
      
      //**CMFD**
      var cmfd_input = dataset.Imgcol_cmfd.filter(filter_date);
      
      if (meth_interp === 'bilinear' || meth_intterp === 'bicubic') {
          cmfd_input = cmfd_input.map(function (img) {
              return img.resample(meth_interp).copyProperties(img, img.propertyNames());
          });
      }
      
      //**LST**
      var LST_daily = dataset.ImgCol_LST.filter(filter_date).map(function (img) {
        return img.subtract(273.15)
                  .copyProperties(img, img.propertyNames()); 
      })
      var cmfd_input_lst = pkg_join.SaveBest(cmfd_input, LST_daily);
      
      // **match CMFD(daily) and modis_input_co2(8 days)**
      var diffFilter = ee.Filter.maxDifference({
        difference: 1000 * 60 * 60 * 24 * 8,
        leftField: 'system:time_start', 
        rightField: 'system:time_start'
      });
      
      var greaterEqFilter = ee.Filter.greaterThanOrEquals({
        leftField: 'system:time_start',
        rightField: 'system:time_start'
      })
      
      var filter = ee.Filter.and(diffFilter, greaterEqFilter)
        
      var pml_input = ee.Join.saveBest('matches', 'measure')
          .apply({
            primary: cmfd_input_lst,
            secondary: modis_input_co2,
            condition: filter
          })
          .map(function(img) { 
            return ee.Image(img).addBands(img.get('matches'))
              .set("matches", null);
          })
      pml_input = ee.ImageCollection(pml_input).sort("system:index")
      
      return ee.ImageCollection(pml_input);
    }
    var INPUTS = PML_INPUTS_daily(year);
    // print(INPUTS)
    
    function PML_period(INPUTS) {
        var len = INPUTS.size();
        /** 2. ImgsRaw: ['Eeq', 'Evp', 'Es_eq', 'Eca', 'Ecr', 'Ei', 'Pi'] */
        var PML_ImgsRaw = INPUTS.map(PML_daily).sort("system:time_start");

        /** 3. Calculate fval_soil, and add Es band */
        var frame = 31; // backward moving average
        var Pi_Es = PML_ImgsRaw.select(['Pi', 'Es_eq']);
        /** movmean_lst(ImgCol, n, win_back = 0, win_forward = 0) */
        var ImgCol_mov = pkg_mov.movmean_lst(Pi_Es, len, frame);
        var fval_soil = ImgCol_mov.map(function (img) {
            return img.expression('b("Pi") / b("Es_eq")').min(1.0).max(0.0)
                .copyProperties(img, pkg_main.global_prop);
        }).select([0], ['fval_soil']);

        /** 4. calculate Es */
        var PML_Imgs_0 = pkg_join.SaveBest(PML_ImgsRaw, fval_soil); //.sort('system:time_start'); 
        var PML_Imgs = PML_Imgs_0.map(function (img) {
            var Es = img.expression('b("Es_eq") * b("fval_soil")').rename('Es').toUint16();
            // var ET = img.expression('b("Ec") + b("Ei") + Es', { Es: Es }).rename('ET');
            return img.addBands(Es); 
        }).select(options.bands.slice(0,5)); 

        return PML_Imgs;
    }
    var PML_Imgs = PML_period(INPUTS);
    
    // Map.addLayer(INPUTS, {}, 'INPUT');
    // Map.addLayer(PML_Imgs, {}, 'PML_Imgs_'.concat(year.toString()));
    return PML_Imgs;
}

pkg_PML.PMLV2_daily = function(){
    var begin_date, ydays;
    
    function contains(xs, x) {
        for (var i = 0; i < xs.length; i++) {
            if (xs[i] === x) return true;
        }
        return false;
    }
    
    for (var year = options.year_begin; year <= options.year_end; year++) {

        var imgcol_PML = PML(year, options.is_PMLV2).aside(print);
        if (options.is_save){
          // define image name for asset
          var dateList = ee.List(imgcol_PML.aggregate_array('system:time_start'))
                          .map(function (date) { return ee.Date(date).format('yyyy-MM-dd'); }).getInfo(); 

          // judge exit imgs in the folder
          if (options.Export.type === "asset") {
              var indexes = ee.ImageCollection(options.Export.folder).aggregate_array('system:index').getInfo();
              // skip already finished
          }
          
          // Batch export 
          var img;
          var n = dateList.length;
          
          for (var i = 0; i < n; i++) {
              var date = dateList[i];
              
              var task = date;
              if (options.Export.type === "asset") {
                if (contains(indexes, task)) continue; // if exist then next
              }
          
              img = imgcol_PML.filterDate(date).first().toUint16().clip(China_10basins); 
              pkg_export.ExportImg(img, task, options.Export);
          }
        }
    }
};

pkg_PML.PMLV2_Yearly = function() {
    var years = pkg_main.seq(options.year_begin, options.year_end);
    
    function func_yearly (year) {
        // year = ee.Number(year);
        var imgcol_PML = PML(year, options.is_PMLV2);
    
        var begin_date = ee.Date.fromYMD(year, 1, 1);
        var task = begin_date.format('YYYY-MM-dd'); //.getInfo();
        var ydays = begin_date.advance(1, 'year').difference(begin_date, 'day');

        var img_year = imgcol_PML.select(options.bands.slice(0, -1)).sum().divide(100) //true values
            .toFloat()
            .set('system:time_start', begin_date.millis())
            .set('system:id', task);
        return img_year;
    }
    
    // enable debug function
    var imgcol_year = years.map(func_yearly);
    
    // print(imgcol_year);
    imgcol_year = ee.ImageCollection(imgcol_year);
    if (options.is_save) pkg_export.ExportImgCol(imgcol_year, options.prefix, options.Export, null);
    return(imgcol_year); 
};

/**
 * [PML_main description]
 *
 * @param {[type]} opt [description]
 * @example

 */
pkg_PML.PML_main = function(opt, verbose) {
    init_dataset();
    pkg_PML.fill_options(opt, verbose);
    
    var imgcol_PML, img_year;

    var debug = false;    
    if (debug) options.is_save = false;

    var imgcol;
    if (options.timescale === "daily") {
        imgcol = pkg_PML.PMLV2_daily();
    } else if (options.timescale === "yearly") {
        imgcol = pkg_PML.PMLV2_Yearly();
        
        if (debug) {
            /** 1. Check the output of PML_V2 **/
            imgcol_PML = PML(2003, options.is_PMLV2);
            print(imgcol_PML, 'imgcol_PML');
            // trend_yearly(imgcol);// yearly trend
            pkg_export.ExportImgCol(imgcol_PML.limit(2), null, options.Export);
        }
    }
    
    // if (verbose) print(imgcol);
    return imgcol;
};

/** 
 * display maps of
 * - trend_GPP
 * - trend_Ec
 * - GPP_2003
 */
 
var vis_et  = { min: 100, max: 1600, palette: pkg_vis.colors.RdYlBu[11] },
    vis_gpp = { min: 100, max: 3500, palette: pkg_vis.colors.RdYlGn[11] };
var vis_slp = { min: -20, max: 20, palette: ["ff0d01", "fafff5", "2aff03"] }; // green to red

var lg_gpp = pkg_vis.grad_legend(vis_gpp, 'GPP', false);
var lg_slp = pkg_vis.grad_legend(vis_slp, 'Trend (gC m-2 y-2)', false); //gC m-2 y-2, kPa y-1
// pkg_vis.add_lgds([lg_gpp, lg_slp]);

function trend_yearly(imgcol_year, show_imgcol){
    show_imgcol = show_imgcol || false;

    // 2. trend
    var img_trend_gpp = pkg_trend.imgcol_trend(imgcol_year, 'GPP', true);
    var img_trend_et = pkg_trend.imgcol_trend(imgcol_year, 'Ec', true);

    Map.addLayer(img_trend_gpp.select('slope'), vis_slp, 'gpp');
    Map.addLayer(img_trend_et.select('slope'), vis_slp, 'Ec');

    var img = imgcol_year.first(); //img_year; //
    Map.addLayer(img.select('GPP'), vis_gpp, 'first_year GPP');
    if (show_imgcol) Map.addLayer(imgcol_year, {}, "imgcol_year");
    print(imgcol_year);
}

pkg_PML.add_ETsum = function(img){
    var ET = img.expression('b("Ec") + b("Ei") + b("Es")').rename("ET");
    return img.addBands(ET);
};

exports = pkg_PML;
/** ------------------------------------------------------------------------ */
var __main__ = true;
if (__main__) {
    // **default is dynamic**
    var opt = {
        year_begin: 2001,
        year_end  : 2001,
        folder    : "XXXX", 
        timescale : "daily", //daily, yearly
        is_dynamic_lc: true, 
        is_save: true
    };
    
    var imgcol_new, imgcol_org;
    
    var imgcol_new = pkg_PML.PML_main(opt, true)//.aside(print).aside(Map.addLayer,{},'PML_new');
    
}
