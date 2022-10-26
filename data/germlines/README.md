
mouse/ is merged from imgt, ogrdb c57bl, and ogrdb balbc

c57bl/ and balbc/ are just from ogrdb

to update mouse germlines, download from: https://ogrdb.airr-community.org/germline_sets/Mouse
  1) cp unaligned fastas to appropriate dir in ogrdb-downloads/
  2) split apart igh files, cp igk/igl files from germlines/mouse/
  3) run script parse-ogrdb.py from partis main dir
  4) after checking that everything went well, replace old mouse/ dirs with new merged-mouse/ dirs
