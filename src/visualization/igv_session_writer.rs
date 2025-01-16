
use log::debug;
use rust_lib_reference_genome::reference_genome::ReferenceGenome;
use simple_error::bail;
use std::collections::BTreeMap;
use std::path::PathBuf;

use crate::data_types::coordinates::Coordinates;
use crate::util::file_io::save_fasta;
use crate::visualization::debug_bam_writer::DebugBamWriter;

// the buffer spacing we put around each side of a consensus sequence
pub const BUFFER_LEN: usize = 1000;
pub const CONTIG_POSTFIX: &str = "custom_contig";

const SESSION_PATH: &str = "custom_igv_session.xml";
const REFERENCE_PATH: &str = "custom_reference.fa";
const REGIONS_PATH: &str = "custom_regions.bed";
const ALIGN_PATH: &str = "custom_alignments.bam";

/// Wrapper stucture that will write out a folder containing everything needed for a custom IGV session
pub struct IgvSessionWriter {
    /// The folder we will save everything to
    session_folder: PathBuf,
    /// Our custom reference genome file
    reference_genome: ReferenceGenome,
    /// Our regions in the reference genome
    regions: Vec<(Coordinates, String)>,
    /// Collection of unmapped records we will map
    unmapped_records: BTreeMap<String, Vec<rust_htslib::bam::Record>>,
    /// If true, adds pre-config on BAM visuals for inspection
    preconfig_bam: bool
}

impl IgvSessionWriter {
    /// Creates a new session writer with no meaningful data. Add data via the `add_custom_region` functionality.
    /// # Arguments
    /// * `session_folder` - the folder we will eventually save everything to
    /// * `preconfig_bam` - if True, this will add some debug rendering options to the BAM XML config
    pub fn new(session_folder: PathBuf, preconfig_bam: bool) -> Self {
        let reference_genome = ReferenceGenome::empty_reference();
        let regions = vec![];
        let unmapped_records = Default::default();
        IgvSessionWriter {
            session_folder,
            reference_genome,
            regions,
            unmapped_records,
            preconfig_bam
        }
    }

    /// Adds a new region to our dataset for mapping
    /// # Arguments
    /// * `region_name` - the contig name to add
    /// * `region_sequence` - the sequence to add, it is recommended to have a `BUFFER_LEN` of Ns at the start / end and in-between anything that should be separated
    /// * `labels` - the labels for this region
    /// * `unmapped_records` - the records that will get mapped to this region
    pub fn add_custom_region(&mut self, region_name: String, region_sequence: &str, region_labels: &[(Coordinates, String)], unmapped_records: Vec<rust_htslib::bam::Record>) -> Result<(), Box<dyn std::error::Error>> {
        // this will throw an error if the region already exists, so no need to check here
        self.reference_genome.add_contig(region_name.clone(), region_sequence)?;

        // add the regions with labels; we need to check the chrom matches on each
        for rl in region_labels.iter() {
            if rl.0.chrom() != region_name {
                bail!("Region {rl:?} is not on correct contig: {region_name}");
            }
        }
        self.regions.extend_from_slice(region_labels);

        // lastly, add the unmapped records with the corresponding contig
        assert!(self.unmapped_records.insert(region_name, unmapped_records).is_none());

        Ok(())
    }

    /// Attempts to save all the data that has been provided to the session writer
    /// # Errors
    /// * if we cannot create the session_folder or any of the sub-files (permissions)
    pub fn write_session(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        // create the folder that captures all our other outputs
        debug!("Creating IGV session folder at {:?}...", self.session_folder);
        std::fs::create_dir_all(&self.session_folder)?;

        // save our reference genome file
        self.save_reference_genome()?;

        // save the regions to a BED
        self.save_regions_bed()?;

        // create our alignments
        self.save_aligned_bam()?;

        // finally save the session file
        self.save_session_file()?;

        Ok(())
    }

    /// Saves the reference genome to a FASTA file
    fn save_reference_genome(&self) -> Result<(), Box<dyn std::error::Error>> {
        let reference_filename = self.session_folder.join(REFERENCE_PATH);
        debug!("Creating custom reference file at {reference_filename:?}");
        
        let mut fasta_map: BTreeMap<String, String> = Default::default();
        for contig_key in self.reference_genome.contig_keys().iter() {
            let value = std::str::from_utf8(self.reference_genome.get_full_chromosome(contig_key))?;
            fasta_map.insert(contig_key.clone(), value.to_string());
        }
        save_fasta(&fasta_map, &reference_filename)?;
        // TODO: ideally, we would bump to latest rust_htslib and do this; however, apparently we get compile issues because 
        //       HiPhase (dep) is on 0.39.5 and I guess that causes problems
        //       likely solution: split off the HiPhase components into a separate crate that can be shared; the components we need do not need htslib
        //       alternate solution: copy the build function in 0.47.0
        // rust_htslib::faidx::build(&reference_filename)?;
        crate::util::file_io::index_fasta(&reference_filename)?;
        Ok(())
    }

    /// This will save the regions of interest to a BED file, which is then overlaid in IGV
    fn save_regions_bed(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        let regions_filename = self.session_folder.join(REGIONS_PATH);
        debug!("Creating custom regions file at {regions_filename:?}");

        let mut bed_writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(regions_filename)?;

        // prior to writing, sort them
        self.regions.sort();

        // chrom start end label
        for (coordinate, label) in self.regions.iter() {
            let start = coordinate.start().to_string();
            let end = coordinate.end().to_string();
            bed_writer.write_record([
                coordinate.chrom(),
                &start,
                &end,
                label
            ])?;
        }

        Ok(())
    }

    /// Creates an aligned BAM file for visualization
    fn save_aligned_bam(&self) -> Result<(), Box<dyn std::error::Error>> {
        let bam_filename = self.session_folder.join(ALIGN_PATH);
        debug!("Creating custom BAM file at {bam_filename:?}");

        // initial the writer with our custom genome
        let mut debug_bam_writer = DebugBamWriter::new(bam_filename, &self.reference_genome)?;

        for (contig_name, unmapped_records) in self.unmapped_records.iter() {
            // map all the reads to it
            let chrom_len = self.reference_genome.get_full_chromosome(contig_name).len();
            let target_region = Coordinates::new(
                contig_name.clone(),
                0_u64,
                chrom_len as u64
            );
            debug_bam_writer.map_records_to_region(unmapped_records, &target_region)?;
        }

        // write all the records we got
        debug_bam_writer.write_all_records()
    }

    /// Create the full session file and saves it
    fn save_session_file(&self) -> Result<(), Box<dyn std::error::Error>> {
        use quick_xml::events::{Event, BytesDecl, BytesEnd, BytesStart};

        let igv_filename = self.session_folder.join(SESSION_PATH);
        debug!("Creating custom IGV session file at {igv_filename:?}");

        let file_handle = std::fs::File::create(igv_filename)?;
        let mut writer = quick_xml::writer::Writer::new_with_indent(file_handle, b' ', 4);
        writer.write_event(Event::Decl(BytesDecl::new("1.0", Some("UTF-8"), Some("no"))))?;

        // write the start of the session block
        let mut session_start = BytesStart::new("Session");
        session_start.push_attribute(("genome", REFERENCE_PATH));
        writer.write_event(Event::Start(session_start))?;

        // this collects all the files as Resources
        writer.write_event(Event::Start(BytesStart::new("Resources")))?;
        
        // add the bam file as a resource
        let mut bam_resource = BytesStart::new("Resource");
        bam_resource.push_attribute(("type", "bam"));
        bam_resource.push_attribute(("path", ALIGN_PATH));
        // bam_resource.push_attribute(("index", format!("{ALIGN_PATH}.bai").as_str())); // for some reason, this actually breaks it
        writer.write_event(Event::Empty(bam_resource))?;

        // add the bed regions as a resource
        let mut bed_resource = BytesStart::new("Resource");
        bed_resource.push_attribute(("type", "bed"));
        bed_resource.push_attribute(("path", REGIONS_PATH));
        writer.write_event(Event::Empty(bed_resource))?;
        writer.write_event(Event::End(BytesEnd::new("Resources")))?;

        // add the first panel, which is the alignment panel
        let mut panel_start = BytesStart::new("Panel");
        panel_start.push_attribute(("name", "Panel0"));
        writer.write_event(Event::Start(panel_start))?;

        let alignment_tracks = [
            vec![
                ("attributeKey", "custom_alignments.bam Coverage"),
                ("autoScale", "true"),
                ("clazz", "org.broad.igv.sam.CoverageTrack"),
                ("id", "custom_alignments.bam_coverage")
            ],
            vec![
                ("attributeKey", "custom_alignments.bam Junctions"),
                ("autoScale", "false"),
                ("clazz", "org.broad.igv.sam.SpliceJunctionTrack"),
                ("id", "custom_alignments.bam_junctions"),
                ("visible", "false")
            ],
            vec![
                ("attributeKey", "custom_alignments.bam"),
                ("clazz", "org.broad.igv.sam.AlignmentTrack"),
                ("id", "custom_alignments.bam")
            ]
        ];
        for track_attributes in alignment_tracks.into_iter() {
            let mut track = BytesStart::new("Track");
            let mut is_bam_track = false;
            for ta in track_attributes.into_iter() {
                track.push_attribute(ta);

                if ta.1 == "custom_alignments.bam" {
                    is_bam_track = true;
                }
            }

            if self.preconfig_bam && is_bam_track {
                // open it
                writer.write_event(Event::Start(track))?;
    
                // add the custom guts: <RenderOptions groupByOption="PHASE" hideSmallIndels="false" quickConsensusMode="false"/>
                let mut render_options = BytesStart::new("RenderOptions");
                render_options.push_attribute(("groupByOption", "PHASE"));
                render_options.push_attribute(("hideSmallIndels", "false"));
                render_options.push_attribute(("quickConsensusMode", "false"));
                writer.write_event(Event::Empty(render_options))?;
    
                // close it
                writer.write_event(Event::End(BytesEnd::new("Track")))?;
    
            } else {
                // just write a basic event
                writer.write_event(Event::Empty(track))?;
            }
        }

        writer.write_event(Event::End(BytesEnd::new("Panel")))?;

        // add the second panel, which has a defined named
        let mut panel_start = BytesStart::new("Panel");
        panel_start.push_attribute(("name", "FeaturePanel"));
        writer.write_event(Event::Start(panel_start))?;

        let feature_tracks = [
            vec![
                ("attributeKey", "Reference sequence"),
                ("clazz", "org.broad.igv.track.SequenceTrack"),
                ("id", "Reference sequence")
            ],
            vec![
                ("attributeKey", "custom_regions.bed"),
                ("clazz", "org.broad.igv.track.FeatureTrack"),
                ("displayMode", "EXPANDED"),
                ("id", "custom_regions.bed")
            ]
        ];
        for track_attributes in feature_tracks.into_iter() {
            let mut track = BytesStart::new("Track");
            for ta in track_attributes.into_iter() {
                track.push_attribute(ta);
            }
            writer.write_event(Event::Empty(track))?;
        }

        writer.write_event(Event::End(BytesEnd::new("Panel")))?;

        // add the layout config
        let mut panel_start = BytesStart::new("PanelLayout");
        panel_start.push_attribute(("dividerFractions", "0.0,0.85"));
        writer.write_event(Event::Empty(panel_start))?;

        // write the end of the session block
        writer.write_event(Event::End(BytesEnd::new("Session")))?;

        Ok(())
    }
}
