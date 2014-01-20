#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/util.hh>
#include <core/fragment/FragmentIO.hh>
#include <protocols/simple_moves/SmoothFragmentMover.hh>
#include <protocols/abinitio/AllResiduesChanged.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/sequence/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/sequence/Sequence.hh>
#include <core/scoring/constraints/util.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/GunnCost.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/SmoothFragmentMover.hh>
#include <protocols/simple_filters/RmsdEvaluator.hh>
#include <protocols/simple_filters/JumpEvaluator.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/WhileMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/simple_filters/RmsdEvaluator.hh>
#include <protocols/simple_filters/JumpEvaluator.hh>
#include <protocols/evaluation/TimeEvaluator.hh>
#include <protocols/evaluation/PCA.hh>
#include <protocols/simple_filters/PoseMetricEvaluator.hh>
#include <protocols/constraints_additional/ConstraintEvaluator.hh>
#include <protocols/evaluation/util.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include "compbio_abinitio.hh"

/* Utility Headers */
#include <core/pose/util.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <numeric/numeric.functions.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/templates.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>
#include <basic/options/keys/filters.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>
#include <basic/options/keys/loopfcst.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/templates.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/after_opts.hh>

//Auto Headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/util.hh>
#include <protocols/moves/MoverStatistics.hh>
#include <ObjexxFCL/format.hh>

//numeric headers
#include <numeric/random/random.hh>

//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL
{
namespace fmt { }
} using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end
using namespace basic;
using namespace core;
using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;

static basic::Tracer tr( "apps.pilot.compbio.compbio_abinitio" );


/* --------------------------register_options---------------------------- */
void
compbio_abinitio::register_options()
{
    Protocol::register_options();
    using namespace basic::options;
    using namespace OptionKeys;
    option.add_relevant( OptionKeys::abinitio::increase_cycles );
    option.add_relevant( OptionKeys::abinitio::smooth_cycles_only );
    option.add_relevant( OptionKeys::abinitio::debug );
    option.add_relevant( OptionKeys::abinitio::skip_convergence_check );
    option.add_relevant( OptionKeys::abinitio::log_frags );
    option.add_relevant( OptionKeys::abinitio::only_stage1 );
    option.add_relevant( OptionKeys::abinitio::end_bias );
    option.add_relevant( OptionKeys::abinitio::symmetry_residue );
    option.add_relevant( OptionKeys::abinitio::vdw_weight_stage1) ;
    option.add_relevant( OptionKeys::abinitio::override_vdw_all_stages );
    option.add_relevant( OptionKeys::abinitio::recover_low_in_stages );
    option.add_relevant( OptionKeys::abinitio::close_chbrk );
    option.add_relevant( OptionKeys::abinitio::bGDT );
    option.add_relevant( OptionKeys::abinitio::rmsd_residues );
    option.add_relevant( OptionKeys::out::file::silent );
    option.add_relevant( OptionKeys::out::nstruct );
}

/* Constructor */
compbio_abinitio::compbio_abinitio( int decoy_num )
    : init_for_input_yet_( false )
{
    temperature_ = 2.0;

    // prepare tag
    std::stringstream oss;
    oss << option[out::output_tag]() << "model_"  << decoy_num;
    tag_ = oss.str();
    tr.Info << "Tag:" << tag_ << std::endl;

    // switching stages on/off
    skip_stage1_ = false;
    skip_stage2_ = false;
    skip_stage3_ = false;
    skip_stage4_ = false;
    skip_stage5_ = true;

    // relevant to stages 2 & 3
    apply_large_frags_ = true;
    short_insert_region_ = false;

    // number of cycles
    num_cycles_stage1_ = num_cycles_stage2_ = num_cycles_stage3_ = 2000;
    num_cycles_stage4_ = 4000;
    num_cycles_stage5_ = 5000;

    // adding evaluator
    using namespace protocols::evaluation;
    evaluator_ = new MetaPoseEvaluator();
}

/* Destructor */
compbio_abinitio::~compbio_abinitio()
{
}

/*---------------------------------run-----------------------------------
 * It first performs initialization, and then folding. At the end, it
 * outputs final 3-D model for protein sequence
 */
void
compbio_abinitio::run()
{
    setup();                      // perform initialization
    init_on_new_input();          // gets wothier in distributed computing

    // print the chosen parameters:
    tr.Info << "num_cycles_stage1 = " <<  num_cycles_stage1_ << std::endl;
    tr.Info << "num_cycles_stage2 = " <<  num_cycles_stage2_ << std::endl;
    tr.Info << "num_cycles_stage3 = " <<  num_cycles_stage3_ << std::endl;
    tr.Info << "num_cycles_stage4 = " <<  num_cycles_stage4_ << std::endl;
    tr.Info << "num_cycles_stage5 = " <<  num_cycles_stage5_ << std::endl;
    tr.Info << "temperature_= " << temperature_ << std::endl;

    fold( *input_pose_ );            // fold the extended chain
    if ( option[out::pdb].user() )
    {
        output_models(); // output final model
    }
}

/*---------------------------------setup-----------------------------------
 * Read 3- and 9-mer fragments, native structure, sequence (either from
 * fasta file or from natve structure, if given), and setting up score
 * functions.
 */
void
compbio_abinitio::setup()
{
    using namespace OptionKeys;
    //using namespace import_pose;
    // Reading 3- and 9-mer fragments
    std::string frag_large_file, frag_small_file;

    if ( option[in::file::frag3].user() ) // smaller 3-mer fragments
    {
        tr.Info << "first, 3-mer fragments ..." << std::endl;
        frag_small_file = option[in::file::frag3]();
        fragset_small_ = fragment::FragmentIO( option[basic::options::OptionKeys::abinitio::number_3mer_frags],
                                               1, option[frags::annotate] ).read_data( frag_small_file );
    }
    if ( option[in::file::frag9].user() ) // larger 9-mer fragments
    {
        tr.Info << "now, 9-mer fragments ..." << std::endl;
        frag_large_file = option[in::file::frag9]();
        fragset_large_ = fragment::FragmentIO( option[basic::options::OptionKeys::abinitio::number_9mer_frags](),
                                               option[frags::nr_large_copies](),
                                               option[frags::annotate]() ).read_data( frag_large_file );
    }

    // Read native pose, if given
    native_pose_ = new pose::Pose;

    if ( option[in::file::native ].user() )
    {
    	using namespace core;
        tr.Info << "Reading native structure ..." << std::endl;
        import_pose::pose_from_pdb( *native_pose_, option[in::file::native]() , false);
        core::pose::set_ss_from_phipsi( *native_pose_ );
        core::util::switch_to_residue_type_set( *native_pose_, chemical::CENTROID );
    }

    // Read target sequence from fasta file or native structure
    if ( option [in::file::fasta].user() )
    {
        sequence_ = core::sequence::read_fasta_file( option[OptionKeys::in::file::fasta]()[1] )[1]->sequence();
        tr.Info << "Read fasta sequence: Size = " << sequence_.size()
                << " Residues = \n"  << sequence_ << std::endl;
    }
    else
    {
        if ( native_pose_ )
        {
            sequence_ = native_pose_->sequence();
            tr.Info << "Sequence taken from native:" << sequence_ << std::endl;
        }
        else
        {
            utility_exit_with_message( "Error: can't read sequence!" );
        }
    }

    tr.Info << "Generate extended pose ...... " << std::endl;
    input_pose_ = new pose::Pose;
    generate_extended_pose( *input_pose_, sequence_ );

    // Setting up score functions
    tr.Info << "Setting up score functions ...... " << std::endl;
    using namespace scoring;
    score_stage1_ = ScoreFunctionFactory::create_score_function( "score0" );
    score_stage2_ = ScoreFunctionFactory::create_score_function( "score1" );
    score_stage3a_ = ScoreFunctionFactory::create_score_function( "score2" );
    score_stage3b_ = ScoreFunctionFactory::create_score_function( "score5" );
    score_stage4_ = ScoreFunctionFactory::create_score_function( "score3" );
    score_stage5_ = ScoreFunctionFactory::create_score_function( "score3" );

    core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn(*score_stage1_);
    core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn(*score_stage2_);
    core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn(*score_stage3a_);
    core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn(*score_stage3b_);
    core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn(*score_stage4_);
    core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn(*score_stage5_);

    // Stages to recover pose
    tr.Info << "Setting update stage recovering vectors ...... " << std::endl;
    recover_low_stages_.clear();
    recover_low_stages_.push_back( STAGE_1 );
    recover_low_stages_.push_back( STAGE_2 );
    recover_low_stages_.push_back( STAGE_3a );
    recover_low_stages_.push_back( STAGE_3b );
    recover_low_stages_.push_back( STAGE_4 );
    recover_low_stages_.push_back( STAGE_5 );

    silent_score_file_ = new io::silent::SilentFileData;
    silent_score_file_->set_filename( option[out::output_tag]() + std::string( option[out::sf]()) );
}

/*---------------------------------fold-------------------------------------
 * Finally start folding the extended protein chain into a three-dimensional
 * structure. In rosetta protocol, the folding process is performed through
 * five stages (involving three kinds of related score functions).
 */
void
compbio_abinitio::fold( core::pose::Pose & pose )
{
    using namespace OptionKeys;
    using namespace moves;

    if ( !init_for_input_yet_ )
        init_on_new_input();

    moves::MonteCarloOP temp_mc_( new moves::MonteCarlo( pose, *score_stage1_, temperature_ ) );
    mc_ = temp_mc_;

    runtime_assert( brute_move_large_ );                      // large
    trial_large_ = new moves::TrialMover( brute_move_large_, mc_ );
    trial_large_->keep_stats_type( all_stats );
    runtime_assert( brute_move_small_ );                     // small
    trial_small_ = new moves::TrialMover( brute_move_small_, mc_ );
    trial_small_->keep_stats_type( all_stats );
    runtime_assert( smooth_move_small_ );                    // smooth small
    smooth_trial_small_ = new moves::TrialMover( smooth_move_small_, mc_ );
    smooth_trial_small_->keep_stats_type( all_stats );

    bool success( true );
    total_trials_ = 0;

    core::io::silent::SilentFileDataOP outsfd( NULL );
    if ( option[out::file::silent].user() )
    {
        outsfd = new core::io::silent::SilentFileData();
    }

    // part 1 ----------------------------------------
    if ( !skip_stage1_ )
    {
        clock_t starttime = clock();

        tr.Info << "\n===================================================\n";
        tr.Info << "Stage 1                                              \n";
        tr.Info << "Folding with score0 for max of " << num_cycles_stage1_ << std::endl;
        if ( option[basic::options::OptionKeys::run::profile] )
            basic::prof_show();

        if ( option[basic::options::OptionKeys::abinitio::debug]() )
            output_debug_structure( pose, "stage0" );

        if ( !get_checkpoints().recover_checkpoint( pose, get_current_tag(),
                "stage_1", false /* fullatom*/, true /*fold tree */) )
        {
            core::scoring::constraints::ConstraintSetOP orig_constraints( NULL );
            orig_constraints = pose.constraint_set()->clone();
            success = fold_stage1( pose );
            if ( tr.Info.visible() )
                current_scorefxn().show( tr, pose );

            recover_low( pose, STAGE_1 );
            mc_->show_counters();
            total_trials_ += mc_->total_trials();
            mc_->reset_counters();
            // restore constraints - this is critical for checkpointing to work
            pose.constraint_set( orig_constraints );
            get_checkpoints().checkpoint( pose, get_current_tag(), "stage_1", true /*fold tree */ );
        } //recover checkpoint

        get_checkpoints().debug( get_current_tag(), "stage_1", current_scorefxn()( pose ) );
        clock_t endtime = clock();

        if ( option [basic::options::OptionKeys::run::profile] )
            basic::prof_show();
        if ( option [basic::options::OptionKeys::abinitio::debug]() )
        {
            tr.Info << "Timeperstep: " << ( double (endtime) - starttime ) / ( CLOCKS_PER_SEC ) << std::endl;
            output_debug_structure( pose, "stage1" );
        }
    } //skipStage1

    if ( !success )
    {
        set_last_move_status( moves::FAIL_RETRY );
        return;
    }

    if ( !skip_stage2_ )
    {
        tr.Info << "\n=====================================================\n";
        tr.Info << " Stage 2                                               \n";
        tr.Info << " Folding with score1 for " << num_cycles_stage2_;
        tr.Info << std::endl;

        clock_t starttime = clock();

        if ( !get_checkpoints().recover_checkpoint( pose, get_current_tag(),
                "stage_2", false /* fullatom */, true /*fold tree */ ) )
        {
            core::scoring::constraints::ConstraintSetOP orig_constraints( NULL );
            orig_constraints = pose.constraint_set()->clone();

            success = fold_stage2( pose );
            recover_low( pose, STAGE_2 );

            if ( tr.visible() )
                current_scorefxn().show( tr, pose );

            mc_->show_counters();
            total_trials_ += mc_->total_trials();
            mc_->reset_counters();

            // restore constraints - this is critical for checkpointing to work
            pose.constraint_set( orig_constraints );
            get_checkpoints().checkpoint( pose, get_current_tag(), "stage_2", true /*fold tree */ );
        }

        get_checkpoints().debug( get_current_tag(), "stage_2", current_scorefxn() (pose) );
        clock_t endtime = clock();

        if ( option [basic::options::OptionKeys::run::profile] )
            basic::prof_show();
        if ( option [basic::options::OptionKeys::abinitio::debug]() )
        {
            output_debug_structure( pose, "stage2" );
            // tr << "Timeperstep: " << ( double (endtime) - starttime ) / ( CLOCKS_PER_SEC ) << std::endl;
        }
    }

    if ( !skip_stage3_ )
    {
        tr.Info << "\n=====================================================\n";
        tr.Info << "   Stage 3                                             \n";
        tr.Info << "   Folding with score2 and score5 for ";
        tr.Info << num_cycles_stage3_ << std::endl;

        clock_t starttime = clock();

        success = fold_stage3( pose );
        recover_low( pose, STAGE_3b );

        if ( tr.Info.visible() )
            current_scorefxn ().show( tr, pose );

        mc_->show_counters();
        total_trials_ += mc_->total_trials();
        mc_->reset_counters();

        clock_t endtime = clock();

        if ( option[basic::options::OptionKeys::run::profile] )
            basic::prof_show ();

        if ( option[ basic::options::OptionKeys::abinitio::debug]() )
        {
            output_debug_structure( pose, "stage3" );
            // tr << "Timeperstep: " << ( double(endtime) - starttime ) / ( CLOCKS_PER_SEC ) << std::endl;
        }
    }

    if ( !skip_stage4_ )
    {
        tr.Info << "\n=====================================================\n";
        tr.Info << "   Stage 4                                             \n";
        tr.Info << "   Folding with score3 for ";
        tr.Info <<  num_cycles_stage4_ << std::endl;

        clock_t starttime = clock();

        success = fold_stage4( pose );
        recover_low( pose, STAGE_4 );

        if ( tr.Info.visible() )
            current_scorefxn().show( tr, pose );

        mc_->show_counters();
        total_trials_ += mc_->total_trials();
        mc_->reset_counters();

        clock_t endtime = clock();

        if ( option[basic::options::OptionKeys::run::profile] )
            basic::prof_show();

        if ( option [basic::options::OptionKeys::abinitio::debug]() )
        {
            output_debug_structure( pose, "stage4" );
            // tr << "Timeperstep: " << ( double (endtime) - starttime ) / ( CLOCKS_PER_SEC ) << std::endl;
        }

        //tr.Info << "\n=====================================================\n";
        //tr.Info << "   Finished Abinitio                                   \n";
        //tr.Info <<  std::endl;
    }

    if ( !skip_stage5_ )
    {
        //tr.Info << "\n=====================================================\n";
        //tr.Info << "   Stage 5                                             \n";
        //tr.Info << "   Folding with score3 for ";
        //tr.Info << num_cycles_stage5_ << std::endl;

        clock_t starttime = clock();

        success = fold_stage5( pose );
        recover_low( pose, STAGE_5 );

        if ( tr.Info.visible() )
            current_scorefxn().show( tr, pose );

        mc_->show_counters();
        total_trials_ += mc_->total_trials();
        mc_->reset_counters();

        clock_t endtime = clock();

        if ( option[basic::options::OptionKeys::run::profile] )
            basic::prof_show();

        if ( option [basic::options::OptionKeys::abinitio::debug]() )
        {
            output_debug_structure( pose, "stage5" );
            // tr << "Timeperstep: " << ( double (endtime) - starttime ) / ( CLOCKS_PER_SEC ) << std::endl;
        }
    }

    tr.Info << "\n=====================================================\n";
    tr.Info << "   Now really finished Abinitio                        \n";
    tr.Info << std::endl;

    get_checkpoints().flush_checkpoints();
    if ( !success )
        set_last_move_status( moves::FAIL_RETRY );

    io::silent::SilentStructOP pss = io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
    core::scoring::ScoreFunction scorefxn = current_scorefxn();
    process_decoy( pose, scorefxn, *pss );

    if ( outsfd )
        outsfd->write_all( option[ out::file::silent ]() );

    if ( silent_score_file_ )
    {
        silent_score_file_->write_silent_struct( *pss, silent_score_file_->filename(), true );
    }
}

void compbio_abinitio::process_decoy( core::pose::Pose & pose, core::scoring::ScoreFunction const & scorefxn,
                                      core::io::silent::SilentStruct & pss ) const
{
    using namespace OptionKeys;
    scorefxn( pose );
    pss.fill_struct( pose, tag_ );

    if ( option[in::file::native ].user() )
    {
        evaluator_->add_evaluation( new protocols::simple_filters::RmsdEvaluator( native_pose_, "_full", option[OptionKeys::abinitio::bGDT]() ) );
    }

    evaluator_->apply( pose, tag_, pss ); // run PoseEvaluators
    if ( option[OptionKeys::jumps::evaluate]() )
    {
        if ( !native_pose_ )
        {
            utility_exit_with_message( "to evaluate jumps, specify native" );
        }
        protocols::evaluation::MetaPoseEvaluatorOP eval_jumps;
        native_pose_->fold_tree( pose.fold_tree() );
        for ( Size nj = 1; nj <= pose.fold_tree().num_jump(); ++nj )
        {
            eval_jumps->add_evaluation( new protocols::simple_filters::JumpEvaluator( *native_pose_, nj ) );
        }
        eval_jumps->apply( pose, tag_, pss );
    }
} // process_decoy

/* ----------------------------output_models------------------------------ */
void
compbio_abinitio::output_models()
{
    std::string pdb_name = tag_ + ".pdb";
    input_pose_->dump_pdb( pdb_name );
}

/* -----------------------------current_scorefxn------------------------------ */
void
compbio_abinitio::current_scorefxn( scoring::ScoreFunction const & scorefxn )
{
    mc_->score_function( scorefxn );
}

/* -----------------------------current_scorefxn------------------------------ */
scoring::ScoreFunction const &
compbio_abinitio::current_scorefxn() const
{
    return mc_->score_function();
}

/* -----------------------------replace_scorefxn------------------------------ */
void
compbio_abinitio::replace_scorefxn( core::pose::Pose & pose, StageID stage )
{
    if ( score_stage1_  && ( stage == STAGE_1 ) )                  // stage 1
        current_scorefxn( *score_stage1_ );
    if ( score_stage2_  && ( stage == STAGE_2 ) )                  // stage 2
        current_scorefxn( *score_stage2_ );
    if ( score_stage3a_ && ( stage == STAGE_3a) )                  // stage 3a
        current_scorefxn( *score_stage3a_ );
    if ( score_stage3b_ && ( stage == STAGE_3b) )                  // stage 3b
        current_scorefxn( *score_stage3b_ );
    if ( score_stage4_  && ( stage == STAGE_4 ) )                  // stage 4
        current_scorefxn( *score_stage4_ );
    if ( score_stage5_  && ( stage == STAGE_5 ) )                  // stage 5
        current_scorefxn( *score_stage5_ );

    Real temperature( temperature_ );
    if ( stage == STAGE_5 )
        temperature = 0.5;
    mc_->set_autotemp( true, temperature );
    mc_->set_temperature( temperature );          // temperature might have
    mc_->reset( pose );                           // changed due to autotemp..

    core::scoring::constraints::add_constraints_from_cmdline_to_pose(pose);
}

/* --------------------------fold_stage1---------------------------- */
bool
compbio_abinitio::fold_stage1( core::pose::Pose & pose )
{
    tr.Debug << "Stage:" << STAGE_1 << "..... " << std::endl;
    replace_scorefxn( pose, STAGE_1 );
    mc_->set_autotemp( false, temperature_ );

    (*score_stage1_)( pose );
    protocols::abinitio::AllResiduesChanged done( pose, brute_move_large_->insert_map(), *movemap_() );
    moves::MoverOP trial( trial_large_ );

    for ( Size j = 1; j <= num_cycles_stage1_; ++j )
    {
        trial->apply( pose );
        if ( done( pose ) )
        {
            tr.Info << "Replaced chain after " << j << " cycles." << std::endl;
            mc_->reset( pose );
            return true;
        }
    }
    tr.Warning << "Warning: extended chain may still remain after ";
    tr.Warning << num_cycles_stage1_ << " cycles!" << std::endl;
    done.show_unmoved( pose, tr.Warning );
    mc_->reset( pose ); // make sure that we keep the final structure

    return true;
}

/* --------------------------fold_stage2---------------------------- */
bool
compbio_abinitio::fold_stage2( core::pose::Pose & pose )
{
    tr.Debug << "Stage:" << STAGE_2 << "..... " << std::endl;
    protocols::abinitio::AllResiduesChanged done( pose, brute_move_large_->insert_map(), *movemap_() );
    moves::TrialMoverOP trials;
    replace_scorefxn( pose, STAGE_2 );
    (*score_stage2_)( pose );

    moves::SequenceMoverOP cycle( new moves::SequenceMover() );
    if ( apply_large_frags_ )
        cycle->add_mover( trial_large_->mover() );
    if ( short_insert_region_ )
        cycle->add_mover( trial_small_->mover() );

    Size nr_cycles = num_cycles_stage2_ / ( short_insert_region_ ? 2 : 1 );
    trials = new moves::TrialMover( cycle, mc_ );
    moves::RepeatMover( trials, nr_cycles ).apply( pose );

    return true;
}

/* --------------------------fold_stage3a---------------------------- */
bool
compbio_abinitio::fold_stage3( core::pose::Pose & pose )
{
    int nloop1 = 1;
    int nloop2 = 10;
    if ( short_insert_region_ )
    {
        nloop1 = 2;
        nloop2 = 5;
    }

    moves::TrialMoverOP trials;
    replace_scorefxn( pose, STAGE_3a );
    (*score_stage3a_)( pose );
    //score for this stage is changed in the do_stage3_cycles explicitly
    if ( option[OptionKeys::templates::change_movemap].user() &&
            option [OptionKeys::templates::change_movemap] == 3 )
    {
        kinematics::MoveMapOP new_mm = new kinematics::MoveMap( *movemap_() );
        new_mm->set_bb( true );
        movemap_ = new_mm;
        if ( smooth_move_small_ )
            smooth_move_small_->set_movemap( movemap_ );
        if ( brute_move_small_ )
            brute_move_small_->set_movemap( movemap_ );
        if ( brute_move_large_ )
            brute_move_large_->set_movemap( movemap_ );
    }

    compbio_convergence_checkOP convergence_checker( NULL );
    if ( !option[OptionKeys::abinitio::skip_convergence_check] )
    {
        convergence_checker = new compbio_convergence_check;
    }
    trials = apply_large_frags_ ? trial_large_ : trial_small_;

    int iteration = 1;
    for ( int lct1 = 1; lct1 <= nloop1; lct1++ )
    {
        if ( lct1 > 1 )
            trials = trial_small_; //only with short_insert_region!
        for ( int lct2 = 1; lct2 <= nloop2; lct2++, iteration++ )
        {
            tr.Debug << "Loop: " << lct1 << "   " << lct2 << std::endl;
            if ( !prepare_loop_in_stage3( pose, iteration) )
                return false;

            if ( !get_checkpoints().recover_checkpoint( pose, get_current_tag(),
                    "stage_3_iter" + ObjexxFCL::string_of( lct1 ) + "_" + string_of( lct2 ),
                    false /*fullatom*/, true /*fold tree*/ ) )
            {
                tr.Debug << "  Score stage3 loop iteration " << lct1 << " " << lct2 << std::endl;
                if ( convergence_checker )
                {
                    moves::TrialMoverOP stage3_trials = trials;
                    convergence_checker->set_trials( stage3_trials );
                    moves::WhileMover( stage3_trials, num_cycles_stage3_, convergence_checker ).apply( pose );
                }
                else
                {
                    moves::RepeatMover( trials, num_cycles_stage3_ ).apply( pose );
                }

                if ( numeric::mod( (int) iteration, 2 ) == 0 || iteration > 7 )
                    recover_low( pose, STAGE_3a );

                recover_low( pose, STAGE_3b );
                get_checkpoints().checkpoint( pose, get_current_tag(),
                                              "stage_3_iter" + string_of( lct1 ) + "_" + string_of( lct2 ), true /*fold tree*/ );
            } //recover_checkpoint
            get_checkpoints().debug( get_current_tag(),
                                     "stage_3_iter" + ObjexxFCL::string_of( lct1 ) + "_" + string_of( lct2 ), current_scorefxn()( pose ) );
        } // loop 2
    } // loop 1
    return true;
}

/* --------------------------fold_stage4---------------------------- */
bool
compbio_abinitio::fold_stage4( core::pose::Pose & pose )
{
    replace_scorefxn( pose, STAGE_4 );
    (*score_stage4_)( pose );
    Size nloop_stage4 = 3;
    moves::TrialMoverOP trials;

    for ( Size kk = 1; kk <= nloop_stage4; ++kk )
    {
        tr.Debug << "prepare:" << std::endl;
        if ( !get_checkpoints().recover_checkpoint( pose, get_current_tag(),
                "stage4_kk_" + ObjexxFCL::string_of( kk ), false /*fullatom*/, true /*fold_tree*/ ) )
        {
            trials = smooth_trial_small_;
            tr.Debug << "start " << num_cycles_stage4_ << " cycles" << std::endl;
            moves::RepeatMover( trials, num_cycles_stage4_ ).apply( pose );
            tr.Debug << "finished" << std::endl;
            recover_low( pose, STAGE_4 );
            get_checkpoints().checkpoint( pose, get_current_tag(),
                                          "stage4_kk_" + ObjexxFCL::string_of( kk ), true /*fold tree*/ );
        }
        get_checkpoints().debug( get_current_tag(),
                                 "stage4_kk_" + ObjexxFCL::string_of( kk ), current_scorefxn()( pose ) );
    }  // loop kk
    return true;
}

/* --------------------------fold_stage5---------------------------- */
bool
compbio_abinitio::fold_stage5( core::pose::Pose & pose )
{
    tr.Debug << "Stage 5: " << STAGE_5 << "..... " << std::endl;

    replace_scorefxn( pose, STAGE_5 );
    (*score_stage5_)( pose );

    Size nmoves = 1;
    moves::TrialMoverOP fa_trials;
    simple_moves::SmallMoverOP fa_small_mover( new simple_moves::SmallMover( movemap_, temperature_, nmoves ) );
    fa_small_mover->angle_max( 'H', 2.0 );
    fa_small_mover->angle_max( 'E', 2.0 );
    fa_small_mover->angle_max( 'L', 5.0 );

    fa_trials = new moves::TrialMover( fa_small_mover, mc_ );
    moves::RepeatMover( fa_trials, num_cycles_stage5_ ).apply( pose );
    mc_->reset( pose );

    return true;
}

/* --------------------------init_on_new_input----------------------------
 * init_on_new_input system allows for initializing these details the first
 * time apply() is called. The job distributor will reinitialize the whole
 * mover when the input changes (a freshly constructed mover, which will
 * re-run this on first apply().
 */
void
compbio_abinitio::init_on_new_input()
{
    protocols::simple_moves::ClassicFragmentMoverOP bms, bml, sms;
    init_for_input_yet_ = true;
    movemap_ = new kinematics::MoveMap;
    movemap_->set_bb( true );
    bms = new protocols::simple_moves::ClassicFragmentMover( fragset_small_, movemap_ );
    bml = new protocols::simple_moves::ClassicFragmentMover( fragset_large_, movemap_ );
    sms = new protocols::simple_moves::SmoothFragmentMover( fragset_small_, movemap_, new protocols::simple_moves::GunnCost );
    bms->set_end_bias( option[OptionKeys::abinitio::end_bias] );
    bml->set_end_bias( option[OptionKeys::abinitio::end_bias] );
    sms->set_end_bias( option[OptionKeys::abinitio::end_bias] );
    brute_move_small_ = bms;
    brute_move_large_ = bml;
    smooth_move_small_ = sms;
}

/* -------------------------prepare_loop_in_stage3----------------------------- */
bool
compbio_abinitio::prepare_loop_in_stage3( core::pose::Pose & pose, Size iteration )
{
    core::Real chbrk_weight_stage_3a = 0;
    core::Real chbrk_weight_stage_3b = 0;
    if ( numeric::mod( (int) iteration, 2 ) == 0 || iteration > 7 )
    {
        Real progress( iteration );
        chbrk_weight_stage_3a = 0.25 * progress;
        tr.Debug << "select score_stage3a..." << std::endl;
        recover_low( pose, STAGE_3a );
        replace_scorefxn( pose, STAGE_3a );
    }
    else
    {
        Real progress( iteration );
        chbrk_weight_stage_3b = 0.05 * progress;
        tr.Debug << "select score_stage3b..." << std::endl;
        recover_low( pose, STAGE_3b );
        replace_scorefxn( pose, STAGE_3b );
    }
    return true;
}

/* -----------------------------contains_stageid----------------------------- */
bool
compbio_abinitio::contains_stageid( utility::vector1< StageID > vec, StageID query )
{
    return find( vec.begin(), vec.end(), query ) != vec.end();
}

/* -----------------------------recover_low---------------------------------- */
void
compbio_abinitio::recover_low( core::pose::Pose & pose, StageID stage )
{
    if ( contains_stageid( recover_low_stages_, stage ) )
    {
        mc_->recover_low( pose );
    }
}

/* -----------------------------generate_extended_pose----------------------------- */
void
compbio_abinitio::generate_extended_pose( core::pose::Pose & pose, std::string const & seq ) const
{
	core::pose::make_pose_from_sequence( pose, seq,
                                       *( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID ) ) );
    for ( Size pos = 1; pos <= pose.total_residue(); pos++ )
    {
        if ( !pose.residue( pos ).is_protein() )
            continue;
        pose.set_phi( pos, -150 );
        pose.set_psi( pos, 150 );
        pose.set_omega( pos, 180 );
    }
}

/* ----------------------------------operator()-------------------------------- */
bool
compbio_convergence_check::operator()( const core::pose::Pose & pose )
{
    if ( !bInit_ )
    {
        bInit_ = true;
        very_old_pose_ = pose;
        return true;
    }
    runtime_assert( trials_ );
    tr.Trace << "TrialCounter in compbio_convergence_check: ";
    tr.Trace << trials_->num_accepts() << std::endl;
    if ( numeric::mod( trials_->num_accepts(), 100 ) != 0 )
        return true;
    if ( (Size) trials_->num_accepts() <= last_move_ )
        return true;
    last_move_ = trials_->num_accepts();
    core::Real converge_rms = core::scoring::CA_rmsd( very_old_pose_, pose );
    very_old_pose_ = pose;
    if ( converge_rms >= 3.0 )
    {
        return true;
    }
    tr.Info << " stop cycles in stage3 due to convergence " << std::endl;

    return false;
}
