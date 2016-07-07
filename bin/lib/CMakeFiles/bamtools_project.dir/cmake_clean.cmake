FILE(REMOVE_RECURSE
  "CMakeFiles/bamtools_project"
  "CMakeFiles/bamtools_project-complete"
  "../externals/bamtools/src/bamtools_project-stamp/bamtools_project-install"
  "../externals/bamtools/src/bamtools_project-stamp/bamtools_project-mkdir"
  "../externals/bamtools/src/bamtools_project-stamp/bamtools_project-download"
  "../externals/bamtools/src/bamtools_project-stamp/bamtools_project-update"
  "../externals/bamtools/src/bamtools_project-stamp/bamtools_project-patch"
  "../externals/bamtools/src/bamtools_project-stamp/bamtools_project-configure"
  "../externals/bamtools/src/bamtools_project-stamp/bamtools_project-build"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/bamtools_project.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
