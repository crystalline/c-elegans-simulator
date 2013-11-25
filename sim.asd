(defsystem wormsim
  :name "wormsim"
  :description "Physical simulation of C.Elegans-like organism"
  :license "MIT"
  :author "Crystalline Emerald"
  :version "0.1"
  :depends-on (:lispbuilder-sdl :lispbuilder-sdl-gfx)
  :serial t
  :components ((:file "sim") (:file "lib3d")))
