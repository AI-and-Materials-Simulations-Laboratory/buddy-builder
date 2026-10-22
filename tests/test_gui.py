"""Test matrix 1 and 16 (headless): widgets drive the backend; reset restores Au/FCC.

The interactive 3D rotation/zoom/pan check needs a browser; see docs/implementation_log.md.
"""

import pytest

pytest.importorskip("ipywidgets")

from nanobuilder.gui import BuilderApp, _uploaded_file  # noqa: E402


@pytest.fixture()
def app():
    return BuilderApp()


def test_launch_state_is_au_fcc_demo(app):
    cfg = app.config()
    assert cfg.material.source == "preset" and cfg.material.preset_id == "Au_fcc"
    assert app.lat_a.value == 4.08
    assert cfg.geometry.shape == "sphere" and cfg.geometry.diameter_nm == 3.0
    assert cfg.ligand.smiles == "SCCO" and cfg.ligand.binding_template == "Au_thiolate"
    assert cfg.coverage.target == 12


def test_generate_enables_downloads_and_edit_invalidates(app):
    assert "No model" in app.downloads.value
    app.generate_btn.click()
    assert app.result.report.status == "ok"
    assert "data:application/zip;base64," in app.downloads.value
    assert "<iframe" in app.viewer.value and "srcdoc=" in app.viewer.value
    app.count.value = 5
    assert "Inputs changed" in app.downloads.value


def test_reset_restores_default_after_upload(app):
    app.example.value = "E_uploaded_mg_cif_slab"
    assert app.config().material.source == "uploaded_file"
    app.reset()
    cfg = app.config()
    assert cfg.material.preset_id == "Au_fcc" and cfg.material.source == "preset"


def test_template_options_follow_material(app):
    assert [k for _, k in app.template.options] == ["Au_thiolate", "chemisorption", "geometric"]
    app.preset.value = "Fe_bcc"
    assert [k for _, k in app.template.options] == ["geometric", "chemisorption"]
    app.preset.value = "Au_fcc"
    assert "Au_thiolate" in [k for _, k in app.template.options]


def test_upload_without_valid_file_does_not_fall_back_to_au(app):
    app.source.value = "uploaded_file"
    app.generate_btn.click()  # goes through the same error handler as a real click
    assert app.result is None
    assert "Upload a crystallographic file" in app.status.value


def test_slab_controls_and_terminations(app):
    app.shape.value = "slab"
    assert [k for _, k in app.region.options] == ["top", "bottom", "both"]
    app.preset.value = "NaCl_rocksalt"
    app.h.value, app.k.value, app.l.value = 1, 1, 1
    assert len(app.termination.options) == 2


def test_visibility_is_per_widget(app):
    assert app.lat_c.layout.display == "none"  # FCC has no free c parameter
    assert app.lat_a.layout.display != "none"
    assert app.diameter.layout.display != "none"
    layouts = [wd.layout for wd in app._config_widgets()]
    assert len({id(lay) for lay in layouts}) == len(layouts)
    app.preset.value = "Mg_hcp"
    assert app.lat_c.layout.display != "none"


def test_file_upload_adapter_handles_both_ipywidgets_versions():
    class Fake:
        pass

    v7, v8 = Fake(), Fake()
    v7.value = {"a.cif": {"content": b"abc", "metadata": {}}}
    v8.value = ({"name": "a.cif", "content": memoryview(b"abc")},)
    assert _uploaded_file(v7) == ("a.cif", b"abc")
    assert _uploaded_file(v8) == ("a.cif", b"abc")


def test_multi_element_crystal_defaults_to_all_anchor_species():
    from nanobuilder.gui import BuilderApp
    app = BuilderApp()
    app.preset.value = "NaCl_rocksalt"
    assert set(app.species_checks) == {"Na", "Cl"}
    assert all(cb.value for cb in app.species_checks.values())
    assert app.species_warning.value == ""
    for cb in app.species_checks.values():
        cb.value = False
    assert "No surface anchor species selected" in app.species_warning.value
    assert app.config().ligand.surface_species == []
    app.species_checks["Na"].value = True
    assert app.species_warning.value == ""
    assert app.config().ligand.surface_species == ["Na"]


def test_library_selection_fills_ligand_inputs():
    from nanobuilder.gui import BuilderApp
    from nanobuilder.ligand_library import LIBRARY
    app = BuilderApp()
    app.library.value = "dodecanethiol"
    cfg = app.config().ligand
    assert cfg.smiles == LIBRARY["dodecanethiol"]["smiles"]
    assert cfg.anchor_atom_index == 0 and cfg.binding_template == "Au_thiolate"
    assert cfg.conformer_mode == "extended" and cfg.library_id == "dodecanethiol"
    app.library.value = "oleylamine"
    cfg = app.config().ligand
    assert cfg.binding_template == "geometric" and abs(cfg.distance_A - 2.4) < 1e-9
    app.smiles.value = "SCCCO"  # a hand edit leaves the library entry
    assert app.library.value == "" and app.config().ligand.library_id is None


def test_thiol_library_entry_on_non_gold_uses_chemisorption():
    from nanobuilder.gui import BuilderApp
    app = BuilderApp()
    app.preset.value = "Ag_fcc"
    app.library.value = "mua"
    assert app.config().ligand.binding_template == "chemisorption"
    assert "chemisorption" in app.library_info.value


def test_gui_rod_controls_reach_the_config():
    from nanobuilder.gui import BuilderApp
    app = BuilderApp()
    app.shape.value = "rod"
    assert [o[1] for o in app.region.options] == ["entire_surface", "side", "tips"]
    assert app.rod_box.layout.display == "" and app.slab_box.layout.display == "none"
    app.rod_len.value, app.rod_caps.value = 12.0, "flat"
    app.rod_u.value, app.rod_v.value, app.rod_w.value = 1, 1, 0
    app.region.value = "tips"
    g = app.config().geometry
    assert (g.shape, g.rod_length_nm, g.rod_caps, g.rod_axis) == ("rod", 12.0, "flat", (1, 1, 0))
    assert "tips" in app.area_info.value


# --- material picker / Materials Project (no network: fake fetch) -------------

def test_typing_a_preset_formula_uses_the_builtin_preset(app):
    app.mat_search.value = "Cu"
    assert app.source.value == "preset" and app.preset.value == "Cu_fcc"


def test_typing_a_library_material_loads_the_preloaded_cif(app):
    app.mat_search.value = "AuCu3"
    assert app.source.value == "materials_project"
    assert app.crystal is not None and set(app.crystal.atoms.get_chemical_symbols()) == {"Au", "Cu"}
    cfg = app.config()
    assert cfg.material.mp_id == app.mp_pick.value and cfg.material.uploaded_bytes
    assert set(app.species_checks) == {"Au", "Cu"} and all(c.value for c in app.species_checks.values())
    app.mat_search.value = "anatase"  # a word of a library name
    assert "anatase" in app.crystal.structure_label
    assert "Space group" in app.cell_info.value and "I4_1/amd" in app.cell_info.value
    app.mat_search.value = "TiO2"  # two library entries: the first loads, the other is named
    assert "rutile" in app.mat_hint.value


def test_unknown_material_is_reported_not_searched_online(app):
    app.mat_search.value = "Unobtainium"
    assert "not in the built-in presets" in app.status.value


def test_browse_menu_has_groups_without_repeated_group_names(app):
    from nanobuilder import material_catalog as catalog
    import json
    assert "bb-groups" in app.mat_menu.value and "Nanorods" in app.mat_menu.value
    assert all(" · " not in label for label in app.mat_search.options)
    label = next(l for l in catalog.labels("Nanorods") if "Gold" in l)
    app.menu_bridge.value = json.dumps({"name": label, "t": 1})
    g = app.config().geometry
    assert app.preset.value == "Au_fcc" and g.shape == "rod" and g.rod_axis == (0, 0, 1)


def test_library_build_and_session_round_trip(app):
    app.mat_search.value = "rutile"
    app.shape.value = "sphere"
    app.diameter.value = 2.0
    app.box_edge.value = 6.0
    app.count.value = 4
    app.generate_btn.click()
    assert app.result.report.ok
    text = (app.package.stage_dir / "session.py").read_text()
    assert "settings(" in text and "build()" in text
    assert (app.package.stage_dir / "certificate_of_build.pdf").read_bytes()[:4] == b"%PDF"
    before = app.config().to_dict()
    app.reset()
    app.console.run(text)
    after = app.config().to_dict()
    assert after["material"]["mp_id"] == before["material"]["mp_id"]
    assert after["geometry"] == before["geometry"] and app.result is not None


def test_theme_switch_sets_root_class_and_rethemes_viewer():
    from nanobuilder.gui import BuilderApp
    app = BuilderApp()
    assert "nb-app" in app.root._dom_classes and "nb-buddy" in app.root._dom_classes
    assert [k for _, k in app.theme.options] == ["buddy", "buddy2", "light", "dark"]
    app.generate_btn.click()
    assert "data-theme=&quot;buddy&quot;" in app.viewer.value
    app.theme.value = "dark"
    assert "nb-dark" in app.root._dom_classes and "nb-buddy" not in app.root._dom_classes
    assert "data-theme=&quot;dark&quot;" in app.viewer.value
    app.theme.value = "light"
    assert "data-theme=&quot;light&quot;" in app.viewer.value


# --- carbon nanomaterials ---------------------------------------------------------


def test_carbon_source_offers_carbon_templates_and_regions(app):
    app.source.value = "carbon"
    assert app.tabs.get_title(2) == "3 Defects & GO" and app.tabs.get_title(1) == "2 Shape"
    assert [k for _, k in app.template.options] == ["C_basal", "C_edge", "C_carboxyl",
                                                    "C_epoxide", "geometric"]
    assert [k for _, k in app.region.options] == ["top", "bottom", "both", "edges",
                                                  "entire_surface"]
    assert app.shape.layout.display == "none" and app.carbon.geometry_box.layout.display == ""
    assert app.carbon.not_carbon.layout.display == "none"
    app.carbon.kind.value = "nanotube"
    assert [k for _, k in app.region.options] == ["entire_surface", "side", "tips", "inner"]
    assert "armchair" in app.carbon.geom_info.value
    app.carbon.kind.value = "carbon_sphere"
    assert [k for _, k in app.region.options] == ["entire_surface"]
    cfg = app.config()
    assert cfg.material.source == "carbon" and cfg.carbon.kind == "carbon_sphere"


def test_carbon_panel_round_trips_its_config(app):
    from nanobuilder.schema import CarbonConfig
    c = CarbonConfig(kind="nanotube", cnt_n=12, cnt_m=4, cnt_walls=2, oxygen_percent=12.5,
                     oxidation_profile="custom", oxidation_pattern="tiles", tile_grid=(2, 4),
                     vacancy_percent=3.5, hole_size="custom", hole_min_atoms=4, hole_max_atoms=9,
                     seed=21)
    app.carbon.apply(c)
    assert app.carbon.config().to_dict() == c.to_dict()


def test_carbon_library_and_covalent_distance(app):
    app.source.value = "carbon"
    app.library.value = "ethanolamine"
    cfg = app.config().ligand
    assert cfg.binding_template == "C_epoxide" and abs(cfg.distance_A - 1.47) < 1e-9
    app.library.value = "phenyl"
    assert app.config().ligand.binding_template == "C_basal"
    assert abs(app.distance.value - 1.53) < 1e-9
    keys = [k for _, k in app.library.options]
    assert "dodecanethiol" not in keys and "cysteamine" in keys  # only ligands that suit carbon
    app.library.value = "cysteamine"  # metal ligand with a carbon route (amide)
    assert app.template.value == "C_carboxyl"
    app.template.value = "C_edge"
    assert abs(app.distance.value - 1.77) < 1e-9  # aryl C-S


def test_typed_carbon_word_opens_the_carbon_entry(app):
    app.mat_search.value = "GO"
    assert app.source.value == "carbon"
    c = app.config().carbon
    assert c.kind == "graphene" and c.oxygen_percent == 30.0 and c.hydroxyl_epoxy_ratio == "2:1"
    app.mat_search.value = "C60"
    assert app.config().carbon.kind == "carbon_sphere"


def test_carbon_generate_preview_and_downloads(app):
    app.source.value = "carbon"
    app.carbon.length.value = app.carbon.width.value = 3.0
    app.carbon.oxygen.value = 20.0
    app.library.value = "ethanolamine"
    app.count.value, app.region.value = 3, "both"
    app.preview_sub_btn.click()
    assert "Oxygen content" in app.status.value and "<iframe" in app.viewer.value
    app.generate_btn.click()
    assert app.result is not None and app.result.report.ok
    assert app.result.assembly.achieved_ligands == 3
    assert "Oxygen content (100 O/C)" in app.results.value
    assert "data:chemical/x-mol2;base64," in app.downloads.value
    app.example.value = "H_cnt_aryl"
    app.generate_btn.click()
    assert app.result.assembly.achieved_ligands == 10
    app.reset()
    assert app.config().material.source == "preset"


# --- Buddy Builder: layout, defaults, channels, boxes, edits, console -------------------


def test_buddy_defaults_basic_mode_and_fibonacci(app):
    header = app.root.children[2].children[0].value
    assert "Buddy Builder" in header and "data:image/png;base64," in header
    assert app.config().coverage.distribution == "fibonacci"  # new builds only
    assert app.mode.value == "basic"
    assert all(box.layout.display == "none" for box in app._advanced_boxes)
    app.mode.value = "advanced"
    assert all(box.layout.display != "none" for box in app._advanced_boxes)
    app.example.value = "B_au_sphere_12_ligands"  # older configs keep their own order
    assert app.config().coverage.distribution == "distributed_random"


def test_material_change_updates_status_and_ligands(app):
    app.source.value = "carbon"
    assert "Material: Graphene" in app.status.value
    assert app.library.value == "carboxyphenyl" and app.template.value == "C_basal"
    assert all(k == "" or k not in ("dodecanethiol", "mercaptoethanol")
               for _, k in app.library.options)
    app.source.value = "preset"
    assert app.library.value == "mercaptoethanol" and app.template.value == "Au_thiolate"
    assert "Material: Gold (Au)" in app.status.value


def test_channel_controls_and_build(app):
    app.shape.value = "channel"
    assert app.channel_box.layout.display != "none" and app.slab_box.layout.display != "none"
    assert [k for _, k in app.region.options][0] == "inner"
    app.channel_type.value = "slit"
    app.channel_count.value = 2
    app.slab_len.value = app.slab_wid.value = 2.0
    app.layers.value = 3
    app.count.value = 8
    app.generate_btn.click()
    assert app.result is not None and app.result.report.ok
    assert app.result.assembly.achieved_ligands == 8
    assert app.config().geometry.channel_count == 2


def test_box_controls_reach_config(app):
    assert app.config().geometry.box_lengths_nm == (8.0, 8.0, 8.0)  # new builds: 8 nm cube
    app.shape.value = "slab"
    assert app.slab_len.layout.display == "none"  # the box sets the slab size
    app.box_shape.value = "rectangular"
    app.box_x.value, app.box_y.value, app.box_z.value = 5.0, 4.0, 6.0
    g = app.config().geometry
    assert g.box_shape == "rectangular" and g.box_lengths_nm == (5.0, 4.0, 6.0)
    app.generate_btn.click()
    cell = app.result.assembly.cell_A
    assert abs(cell[2][2] - 60.0) < 1e-6 and abs(cell[0][0] - 50.0) < 3.0
    app.box_shape.value = "auto"
    assert app.config().geometry.box_lengths_nm is None and app.slab_len.layout.display != "none"


def test_sad_buddy_on_failed_build(app):
    app.box_edge.value = 2.0  # a 3 nm sphere does not fit a 2 nm cube
    app.generate_btn.click()
    assert app.result is None and "enlarge the box" in app.status.value
    assert "bb-sad" in app.popup.value and "Oh no! Buddy could not build that." in app.popup.value
    assert "How to fix" in app.popup.value and "How to fix" in app.status.value
    assert "make the box larger" in app.popup.value


def test_save_as_names_the_downloads(app):
    app.save_as.value = "my gold ball!"
    app.generate_btn.click()
    assert app.package.zip_name == "my_gold_ball.zip"
    assert 'download="my_gold_ball_certificate_of_build.pdf"' in app.downloads.value
    assert 'download="my_gold_ball_session.py"' in app.downloads.value
    app.save_as.value = ""
    app.generate_btn.click()
    from nanobuilder.reporting import ENDEARMENTS
    assert any(app.package.zip_name.startswith(e + "_Buddy_") for e in ENDEARMENTS)


def test_viewer_edits_and_undo(app):
    import json
    app.generate_btn.click()
    assert '&quot;editable&quot;:true' in app.viewer.value
    app.bridge.value = json.dumps({"action": "remove", "ligand": 2, "t": 1})
    assert app.result.assembly.achieved_ligands == 11
    assert "flagged in README" in app.edit_info.value and not app.undo_btn.disabled
    readme = (app.package.stage_dir / "README.txt").read_text()
    assert "MANUAL EDITS" in readme
    app.undo_btn.click()
    assert app.result.assembly.achieved_ligands == 12


def test_console_commands_drive_the_gui(app):
    app.console.run('material(preset="Au_fcc"); shape("sphere", diameter=2.5)')
    assert app.diameter.value == 2.5
    app.console.run('ligand("dodecanethiol"); coverage(count=6, order="even")')
    assert app.config().ligand.library_id == "dodecanethiol"
    assert app.config().coverage.distribution == "evenly_spaced"
    app.console.run("build()")
    assert app.result.assembly.achieved_ligands == 6
    app.console.run('select(element="S")')
    assert len(app.console.selection) == 6
    app.console.run("measure(0, 1)")
    assert "distance 0-1" in app.console_log.value
    app.console.run("__import__('os')")  # not a command: refused, nothing runs
    assert "Unknown command" in app.console_log.value
    app.console.run("x = 1")
    assert "Only commands" in app.console_log.value
