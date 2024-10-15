use crate::constants::*;
use crate::structs::variant::Variant;
use crate::structs::variant_call::VariantCall;
use crate::structs::variants_list::VariantsList;


/// Verifies that `VariantsList::merge` function merges two of the same insertion variants.
///
/// # Scenario
/// One variants list has two variants `chr1:100:SNV:T>A` and `chr1:1000:INS:AAA`.
/// Another variants list has one variant `chr1:1000:INS:AAA`.
///
/// # Expected
/// The merged variants list should have 2 variants `chr1:100:SNV:T>A` and `chr1:1000:INS:AAA`.
#[test]
fn merge_variants_lists_1() {
    let mut variants_list_1: VariantsList = VariantsList::new();
    let mut variants_list_2: VariantsList = VariantsList::new();
    let mut variant_1: Variant = Variant::new("variant_1".to_string());
    let mut variant_2: Variant = Variant::new("variant_2".to_string());
    let mut variant_3: Variant = Variant::new("variant_3".to_string());
    let variant_call_1: VariantCall = VariantCall::new(
        "variant_call_1".to_string(),
        "sample".to_string(),
        "chr1".to_string(),
        100,
        "chr1".to_string(),
        100,
        SINGLE_NUCLEOTIDE_VARIANT.to_string(),
        "T".to_string(),
        "A".to_string(),
        "source".to_string(),
        "".to_string(),
        "".to_string(),
        "DNA".to_string(),
        "method".to_string(),
        "platform".to_string(),
        "PASS".to_string(),
        60.0,
        "yes".to_string(),
        "".to_string(),
        1,
        3,
        3,
        6,
        0.5,
        1000,
        60.0,
        60.0
    );
    let mut variant_call_2: VariantCall = variant_call_1.clone();
    variant_call_2.id = "variant_call_2".to_string();
    variant_call_2.variant_type = INSERTION.to_string();
    variant_call_2.reference_allele = "T".to_string();
    variant_call_2.alternate_allele = "TAAA".to_string();
    variant_call_2.variant_size = 3;
    variant_call_2.position_1 = 1000;
    variant_call_2.position_2 = 1000;
    let mut variant_call_3: VariantCall = variant_call_2.clone();
    variant_call_3.id = "variant_call_3".to_string();

    variant_1.add_variant_call(variant_call_1);
    variant_2.add_variant_call(variant_call_2);
    variant_3.add_variant_call(variant_call_3);
    variants_list_1.add_variant(variant_1);
    variants_list_1.add_variant(variant_2);
    variants_list_2.add_variant(variant_3);

    let mut variants_lists: Vec<&VariantsList> = Vec::new();
    variants_lists.push(&variants_list_1);
    variants_lists.push(&variants_list_2);
    let merged_variants_list: VariantsList = VariantsList::merge(
        &variants_lists,
        2,
        10,
        true,
        true,
        0.5,
        0.5,
        &VARIANT_TYPES_MAP
    );
    assert_eq!(merged_variants_list.variants.len(), 2);
}

/// Verifies that `VariantsList::merge` function does not merge two insertions in the same
/// position with significantly different sizes.
///
/// # Scenario
/// One variants list has two variants `chr1:100:SNV:T>A` and `chr1:1000:INS:TAAAATCACGATCATCATCGACTACTACGTG`.
/// Another variants list has one variant `chr1:1000:INS:AAA`.
///
/// # Expected
/// The merged variants list should have 3 variants `chr1:100:SNV:T>A`, `chr1:1000:INS:AAA`, and
/// `chr1:1000:INS:TAAAATCACGATCATCATCGACTACTACGTG`.
#[test]
fn merge_variants_lists_2() {
    let mut variants_list_1: VariantsList = VariantsList::new();
    let mut variants_list_2: VariantsList = VariantsList::new();
    let mut variant_1: Variant = Variant::new("variant_1".to_string());
    let mut variant_2: Variant = Variant::new("variant_2".to_string());
    let mut variant_3: Variant = Variant::new("variant_3".to_string());
    let variant_call_1: VariantCall = VariantCall::new(
        "variant_call_1".to_string(),
        "sample".to_string(),
        "chr1".to_string(),
        100,
        "chr1".to_string(),
        100,
        SINGLE_NUCLEOTIDE_VARIANT.to_string(),
        "T".to_string(),
        "A".to_string(),
        "source".to_string(),
        "".to_string(),
        "".to_string(),
        "DNA".to_string(),
        "method".to_string(),
        "platform".to_string(),
        "PASS".to_string(),
        60.0,
        "yes".to_string(),
        "".to_string(),
        1,
        3,
        3,
        6,
        0.5,
        1000,
        60.0,
        60.0
    );
    let mut variant_call_2: VariantCall = variant_call_1.clone();
    variant_call_2.id = "variant_call_2".to_string();
    variant_call_2.variant_type = INSERTION.to_string();
    variant_call_2.reference_allele = "T".to_string();
    variant_call_2.alternate_allele = "TAAAATCACGATCATCATCGACTACTACGTG".to_string();
    variant_call_2.variant_size = 30;
    variant_call_2.position_1 = 1000;
    variant_call_2.position_2 = 1000;
    let mut variant_call_3: VariantCall = variant_call_2.clone();
    variant_call_3.id = "variant_call_3".to_string();
    variant_call_3.alternate_allele = "TAAA".to_string();
    variant_call_3.variant_size = 3;

    variant_1.add_variant_call(variant_call_1);
    variant_2.add_variant_call(variant_call_2);
    variant_3.add_variant_call(variant_call_3);
    variants_list_1.add_variant(variant_1);
    variants_list_1.add_variant(variant_2);
    variants_list_2.add_variant(variant_3);

    let mut variants_lists: Vec<&VariantsList> = Vec::new();
    variants_lists.push(&variants_list_1);
    variants_lists.push(&variants_list_2);
    let merged_variants_list: VariantsList = VariantsList::merge(
        &variants_lists,
        2,
        10,
        true,
        true,
        0.5,
        0.5,
        &VARIANT_TYPES_MAP
    );
    assert_eq!(merged_variants_list.variants.len(), 3);
}

/// Verifies that `VariantsList::merge` function merges two deletions in adjacent positions
/// with similar sizes.
///
/// # Scenario
/// One variants list has two variants `chr1:100:SNV:T>A` and `chr1:1001-1030:DEL`.
/// Another variants list has one variant `chr1:1001-1031:DEL`.
///
/// # Expected
/// The merged variants list should have 2 variants `chr1:100:SNV:T>A`,
/// `[chr1:1001-1030:DEL, chr1:1001-1031:DEL]`.
#[test]
fn merge_variants_lists_3() {
    let mut variants_list_1: VariantsList = VariantsList::new();
    let mut variants_list_2: VariantsList = VariantsList::new();
    let mut variant_1: Variant = Variant::new("variant_1".to_string());
    let mut variant_2: Variant = Variant::new("variant_2".to_string());
    let mut variant_3: Variant = Variant::new("variant_3".to_string());
    let variant_call_1: VariantCall = VariantCall::new(
        "variant_call_1".to_string(),
        "sample".to_string(),
        "chr1".to_string(),
        100,
        "chr1".to_string(),
        100,
        SINGLE_NUCLEOTIDE_VARIANT.to_string(),
        "T".to_string(),
        "A".to_string(),
        "source".to_string(),
        "".to_string(),
        "".to_string(),
        "DNA".to_string(),
        "method".to_string(),
        "platform".to_string(),
        "PASS".to_string(),
        60.0,
        "yes".to_string(),
        "".to_string(),
        1,
        3,
        3,
        6,
        0.5,
        1000,
        60.0,
        60.0
    );
    let mut variant_call_2: VariantCall = variant_call_1.clone();
    variant_call_2.id = "variant_call_2".to_string();
    variant_call_2.variant_type = DELETION.to_string();
    variant_call_2.reference_allele = "TAAAATCACGATCATCATCGACTACTACGT".to_string();
    variant_call_2.alternate_allele = "".to_string();
    variant_call_2.variant_size = 30;
    variant_call_2.position_1 = 1001;
    variant_call_2.position_2 = 1030;
    let mut variant_call_3: VariantCall = variant_call_2.clone();
    variant_call_3.id = "variant_call_3".to_string();
    variant_call_2.reference_allele = "TAAAATCACGATCATCATCGACTACTACGTG".to_string();
    variant_call_2.alternate_allele = "".to_string();
    variant_call_2.variant_size = 31;
    variant_call_2.position_1 = 1001;
    variant_call_2.position_2 = 1031;

    variant_1.add_variant_call(variant_call_1);
    variant_2.add_variant_call(variant_call_2);
    variant_3.add_variant_call(variant_call_3);
    variants_list_1.add_variant(variant_1);
    variants_list_1.add_variant(variant_2);
    variants_list_2.add_variant(variant_3);

    let mut variants_lists: Vec<&VariantsList> = Vec::new();
    variants_lists.push(&variants_list_1);
    variants_lists.push(&variants_list_2);
    let merged_variants_list: VariantsList = VariantsList::merge(
        &variants_lists,
        2,
        10,
        true,
        true,
        0.5,
        0.5,
        &VARIANT_TYPES_MAP
    );
    assert_eq!(merged_variants_list.variants.len(), 2);
}

/// Verifies that `VariantsList::merge` function does not merge deletions in adjacent positions
/// with different sizes.
///
/// # Scenario
/// One variants list has two variants `chr1:100:SNV:T>A` and `chr1:1001-1030:DEL`.
/// Another variants list has one variant `chr1:1001-1031:DEL`.
///
/// # Expected
/// The merged variants list should have 3 variants `chr1:100:SNV:T>A`,
/// `chr1:1001-1030:DEL`, and `chr1:1001-1008:DEL`.
#[test]
fn merge_variants_lists_4() {
    let mut variants_list_1: VariantsList = VariantsList::new();
    let mut variants_list_2: VariantsList = VariantsList::new();
    let mut variant_1: Variant = Variant::new("variant_1".to_string());
    let mut variant_2: Variant = Variant::new("variant_2".to_string());
    let mut variant_3: Variant = Variant::new("variant_3".to_string());
    let variant_call_1: VariantCall = VariantCall::new(
        "variant_call_1".to_string(),
        "sample".to_string(),
        "chr1".to_string(),
        100,
        "chr1".to_string(),
        100,
        SINGLE_NUCLEOTIDE_VARIANT.to_string(),
        "T".to_string(),
        "A".to_string(),
        "source".to_string(),
        "".to_string(),
        "".to_string(),
        "DNA".to_string(),
        "method".to_string(),
        "platform".to_string(),
        "PASS".to_string(),
        60.0,
        "yes".to_string(),
        "".to_string(),
        1,
        3,
        3,
        6,
        0.5,
        1000,
        60.0,
        60.0
    );
    let mut variant_call_2: VariantCall = variant_call_1.clone();
    variant_call_2.id = "variant_call_2".to_string();
    variant_call_2.variant_type = DELETION.to_string();
    variant_call_2.reference_allele = "TAAAATCACGATCATCATCGACTACTACGT".to_string();
    variant_call_2.alternate_allele = "".to_string();
    variant_call_2.variant_size = 30;
    variant_call_2.position_1 = 1001;
    variant_call_2.position_2 = 1030;
    let mut variant_call_3: VariantCall = variant_call_2.clone();
    variant_call_3.id = "variant_call_3".to_string();
    variant_call_2.reference_allele = "TAAAATCA".to_string();
    variant_call_2.alternate_allele = "".to_string();
    variant_call_2.variant_size = 8;
    variant_call_2.position_1 = 1001;
    variant_call_2.position_2 = 1008;

    variant_1.add_variant_call(variant_call_1);
    variant_2.add_variant_call(variant_call_2);
    variant_3.add_variant_call(variant_call_3);
    variants_list_1.add_variant(variant_1);
    variants_list_1.add_variant(variant_2);
    variants_list_2.add_variant(variant_3);

    let mut variants_lists: Vec<&VariantsList> = Vec::new();
    variants_lists.push(&variants_list_1);
    variants_lists.push(&variants_list_2);
    let merged_variants_list: VariantsList = VariantsList::merge(
        &variants_lists,
        2,
        10,
        true,
        true,
        0.5,
        0.5,
        &VARIANT_TYPES_MAP
    );
    assert_eq!(merged_variants_list.variants.len(), 3);
}

/// Verifies that `VariantsList::merge` function merges translocations with adjacent
/// breakpoint positions.
///
/// # Scenario
/// One variants list has two variants `chr1:100:SNV:T>A` and `chr1:50000-chr2:50000:BND`.
/// Another variants list has one variant `chr1:50500-chr2:49500:BND`.
///
/// # Expected
/// The merged variants list should have 2 variants `chr1:100:SNV:T>A` and
/// `[chr1:50000-chr2:50000:BND, chr1:50500-chr2:49500:BND]`.
#[test]
fn merge_variants_lists_5() {
    let mut variants_list_1: VariantsList = VariantsList::new();
    let mut variants_list_2: VariantsList = VariantsList::new();
    let mut variant_1: Variant = Variant::new("variant_1".to_string());
    let mut variant_2: Variant = Variant::new("variant_2".to_string());
    let mut variant_3: Variant = Variant::new("variant_3".to_string());
    let variant_call_1: VariantCall = VariantCall::new(
        "variant_call_1".to_string(),
        "sample".to_string(),
        "chr1".to_string(),
        100,
        "chr1".to_string(),
        100,
        SINGLE_NUCLEOTIDE_VARIANT.to_string(),
        "T".to_string(),
        "A".to_string(),
        "source".to_string(),
        "".to_string(),
        "".to_string(),
        "DNA".to_string(),
        "method".to_string(),
        "platform".to_string(),
        "PASS".to_string(),
        60.0,
        "yes".to_string(),
        "".to_string(),
        1,
        3,
        3,
        6,
        0.5,
        1000,
        60.0,
        60.0
    );
    let mut variant_call_2: VariantCall = variant_call_1.clone();
    variant_call_2.id = "variant_call_2".to_string();
    variant_call_2.chromosome_1 = "chr1".to_string();
    variant_call_2.chromosome_2 = "chr2".to_string();
    variant_call_2.variant_type = BREAKPOINT.to_string();
    variant_call_2.reference_allele = "".to_string();
    variant_call_2.alternate_allele = "".to_string();
    variant_call_2.variant_size = -1;
    variant_call_2.position_1 = 50000;
    variant_call_2.position_2 = 50000;
    let mut variant_call_3: VariantCall = variant_call_2.clone();
    variant_call_3.id = "variant_call_3".to_string();
    variant_call_2.position_1 = 50500;
    variant_call_2.position_2 = 49500;

    variant_1.add_variant_call(variant_call_1);
    variant_2.add_variant_call(variant_call_2);
    variant_3.add_variant_call(variant_call_3);
    variants_list_1.add_variant(variant_1);
    variants_list_1.add_variant(variant_2);
    variants_list_2.add_variant(variant_3);

    let mut variants_lists: Vec<&VariantsList> = Vec::new();
    variants_lists.push(&variants_list_1);
    variants_lists.push(&variants_list_2);
    let merged_variants_list: VariantsList = VariantsList::merge(
        &variants_lists,
        2,
        1000,
        true,
        true,
        0.5,
        0.5,
        &VARIANT_TYPES_MAP
    );
    assert_eq!(merged_variants_list.variants.len(), 2);
}

/// Verifies that `VariantsList::merge` function does not merge translocations that share
/// only one of their breakpoint positions.
///
/// # Scenario
/// One variants list has two variants `chr1:100:SNV:T>A` and `chr1:50000-chr2:50000:BND`.
/// Another variants list has one variant `chr1:50500-chr3:80000:BND`.
///
/// # Expected
/// The merged variants list should have 3 variants `chr1:100:SNV:T>A` and
/// `chr1:50000-chr2:50000:BND`, and `chr1:50500-chr3:80000:BND`.
#[test]
fn merge_variants_lists_6() {
    let mut variants_list_1: VariantsList = VariantsList::new();
    let mut variants_list_2: VariantsList = VariantsList::new();
    let mut variant_1: Variant = Variant::new("variant_1".to_string());
    let mut variant_2: Variant = Variant::new("variant_2".to_string());
    let mut variant_3: Variant = Variant::new("variant_3".to_string());
    let variant_call_1: VariantCall = VariantCall::new(
        "variant_call_1".to_string(),
        "sample".to_string(),
        "chr1".to_string(),
        100,
        "chr1".to_string(),
        100,
        SINGLE_NUCLEOTIDE_VARIANT.to_string(),
        "T".to_string(),
        "A".to_string(),
        "source".to_string(),
        "".to_string(),
        "".to_string(),
        "DNA".to_string(),
        "method".to_string(),
        "platform".to_string(),
        "PASS".to_string(),
        60.0,
        "yes".to_string(),
        "".to_string(),
        1,
        3,
        3,
        6,
        0.5,
        1000,
        60.0,
        60.0
    );
    let mut variant_call_2: VariantCall = variant_call_1.clone();
    variant_call_2.id = "variant_call_2".to_string();
    variant_call_2.chromosome_1 = "chr1".to_string();
    variant_call_2.chromosome_2 = "chr2".to_string();
    variant_call_2.variant_type = BREAKPOINT.to_string();
    variant_call_2.reference_allele = "".to_string();
    variant_call_2.alternate_allele = "".to_string();
    variant_call_2.variant_size = -1;
    variant_call_2.position_1 = 50000;
    variant_call_2.position_2 = 50000;
    let mut variant_call_3: VariantCall = variant_call_2.clone();
    variant_call_3.id = "variant_call_3".to_string();
    variant_call_2.chromosome_2 = "chr3".to_string();
    variant_call_2.position_1 = 50500;
    variant_call_2.position_2 = 80000;

    variant_1.add_variant_call(variant_call_1);
    variant_2.add_variant_call(variant_call_2);
    variant_3.add_variant_call(variant_call_3);
    variants_list_1.add_variant(variant_1);
    variants_list_1.add_variant(variant_2);
    variants_list_2.add_variant(variant_3);

    let mut variants_lists: Vec<&VariantsList> = Vec::new();
    variants_lists.push(&variants_list_1);
    variants_lists.push(&variants_list_2);
    let merged_variants_list: VariantsList = VariantsList::merge(
        &variants_lists,
        2,
        1000,
        false,
        true,
        0.5,
        0.5,
        &VARIANT_TYPES_MAP
    );
    assert_eq!(merged_variants_list.variants.len(), 2);
}

/// Verifies that `VariantsList::merge` function merges two of three  variants that are all
/// close to one another but two of the three are of similar size.
///
/// # Scenario
/// One variants list has two variants `chr1:100:INS:A` and `chr1:150:INS:GATACGATAA`.
/// Another variants list has one variant `chr1:145:INS:GATACGATAA`.
///
/// # Expected
/// The merged variants list should have 2 variant2 `chr1:100:INS:A` and
/// `[chr1:110:INS:TGATACGATAA, chr1:145:INS:GATACGATAA]`.
#[test]
fn merge_variants_lists_7() {
    let mut variants_list_1: VariantsList = VariantsList::new();
    let mut variants_list_2: VariantsList = VariantsList::new();
    let mut variant_1: Variant = Variant::new("variant_1".to_string());
    let mut variant_2: Variant = Variant::new("variant_2".to_string());
    let mut variant_3: Variant = Variant::new("variant_3".to_string());
    let variant_call_1: VariantCall = VariantCall::new(
        "variant_call_1".to_string(),
        "sample".to_string(),
        "chr1".to_string(),
        100,
        "chr1".to_string(),
        100,
        INSERTION.to_string(),
        "T".to_string(),
        "TA".to_string(),
        "source".to_string(),
        "".to_string(),
        "".to_string(),
        "DNA".to_string(),
        "method".to_string(),
        "platform".to_string(),
        "PASS".to_string(),
        60.0,
        "yes".to_string(),
        "".to_string(),
        1,
        3,
        3,
        6,
        0.5,
        1000,
        60.0,
        60.0
    );
    let mut variant_call_2: VariantCall = variant_call_1.clone();
    variant_call_2.id = "variant_call_2".to_string();
    variant_call_2.reference_allele = "T".to_string();
    variant_call_2.alternate_allele = "TGATACGATAA".to_string();
    variant_call_2.variant_size = 10;
    variant_call_2.position_1 = 150;
    variant_call_2.position_2 = 150;
    let mut variant_call_3: VariantCall = variant_call_2.clone();
    variant_call_3.id = "variant_call_3".to_string();
    variant_call_2.position_1 = 145;
    variant_call_2.position_2 = 145;
    variant_call_2.reference_allele = "T".to_string();
    variant_call_2.alternate_allele = "TGATACGATAA".to_string();

    variant_1.add_variant_call(variant_call_1);
    variant_2.add_variant_call(variant_call_2);
    variant_3.add_variant_call(variant_call_3);
    variants_list_1.add_variant(variant_1);
    variants_list_1.add_variant(variant_2);
    variants_list_2.add_variant(variant_3);

    let mut variants_lists: Vec<&VariantsList> = Vec::new();
    variants_lists.push(&variants_list_1);
    variants_lists.push(&variants_list_2);
    let merged_variants_list: VariantsList = VariantsList::merge(
        &variants_lists,
        2,
        100,
        true,
        true,
        0.5,
        0.5,
        &VARIANT_TYPES_MAP
    );

    assert_eq!(merged_variants_list.variants.len(), 2);
}