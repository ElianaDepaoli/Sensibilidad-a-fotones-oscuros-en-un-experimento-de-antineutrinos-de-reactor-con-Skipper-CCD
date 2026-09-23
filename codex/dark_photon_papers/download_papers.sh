#!/usr/bin/env bash
set -u

paper_dir="$(cd "$(dirname "$0")" && pwd)"

download_paper() {
  local file_name="$1"
  local arxiv_id="$2"
  if [[ -s "$paper_dir/$file_name" ]]; then
    return 0
  fi
  curl -L --fail --retry 2 --output "$paper_dir/$file_name" "https://arxiv.org/pdf/$arxiv_id"
}

download_paper 01_fabbrichesi_2020_dark_photon_review.pdf 2005.01515
download_paper 02_graham_2021_accelerator_review.pdf 2104.10280
download_paper 03_essig_2013_dark_sectors_snowmass.pdf 1311.0029
download_paper 04_bjorken_2009_fixed_target_theory.pdf 0906.0580
download_paper 05_arias_2012_wispy_cold_dark_matter.pdf 1201.5902
download_paper 06_graham_2016_vector_dm_inflation.pdf 1504.02102
download_paper 07_an_2013_stellar_constraints.pdf 1302.3884
download_paper 08_redondo_2013_solar_constraints.pdf 1305.2920
download_paper 09_chang_2017_sn1987a.pdf 1611.03864
download_paper 10_babar_2014_visible.pdf 1406.2980
download_paper 11_babar_2017_invisible.pdf 1702.03327
download_paper 12_na48_2015_pi0.pdf 1504.00607
download_paper 13_na64_2016_invisible.pdf 1610.02988
download_paper 14_lhcb_2017_dimuon.pdf 1710.02867
download_paper 15_faser_2023_dark_photon.pdf 2308.05587
download_paper 16_dark_srf_2023_lsw.pdf 2301.11512
download_paper 17_dixit_2020_qubit.pdf 2008.12231
download_paper 18_cervantes_2022_srf_haloscope.pdf 2208.03183
download_paper 19_admx_orpheus_2022.pdf 2204.09475
download_paper 20_dosue_rr_2022.pdf 2205.03679
download_paper 21_gigabread_2023.pdf 2310.13891
download_paper 22_ldmx_2018_design.pdf 1808.05219
download_paper 23_darkquest_2022.pdf 2203.08322
download_paper 24_belle2_2022_long_lived.pdf 2202.03452
download_paper 25_jiang_2024_sensor_network.pdf 2305.00890
download_paper 26_nasduck_2026_magnetometry.pdf 2602.22308
download_paper 27_caputo_essig_2026_perspective.pdf 2603.08430
download_paper 28_ldmx_2026_long_lived_sensitivity.pdf 2604.14359
