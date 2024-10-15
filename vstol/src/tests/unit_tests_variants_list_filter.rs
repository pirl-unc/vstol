use serde_json::Value;
use crate::constants::*;
use crate::structs::variant::Variant;
use crate::structs::variant_call::VariantCall;
use crate::structs::variant_filter::VariantFilter;
use crate::structs::variants_list::VariantsList;

/// Verifies that `VariantsList::filter` function filters variants by `filter`.
///
/// # Scenario
/// Filter a variant without PASS filter.
///
/// # Expected
/// The variants list should have 1 variant.
#[test]
fn filter_variants_list_1() {
    let mut variants_list: VariantsList = VariantsList::new();
    let mut variant_1: Variant = Variant::new("variant_1".to_string());
    let mut variant_2: Variant = Variant::new("variant_2".to_string());
    let variant_call_1: VariantCall = VariantCall::new(
        "variant_call_1".to_string(),
        "sample".to_string(),
        "chr1".to_string(),
        100,
        "chr1".to_string(),
        100,
        INSERTION.to_string(),
        "T".to_string(),
        "TAAA".to_string(),
        "source".to_string(),
        "".to_string(),
        "".to_string(),
        "DNA".to_string(),
        "method".to_string(),
        "platform".to_string(),
        "germline".to_string(),
        60.0,
        "yes".to_string(),
        "".to_string(),
        3,
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
    variant_call_2.position_1 = 1000;
    variant_call_2.position_2 = 1000;
    variant_call_2.filter = "PASS".to_string();
    variant_1.add_variant_call(variant_call_1);
    variant_2.add_variant_call(variant_call_2);
    variants_list.add_variant(variant_1);
    variants_list.add_variant(variant_2);

    let mut variant_filters: Vec<VariantFilter> = Vec::new();
    let variant_filter: VariantFilter = VariantFilter::new(
        QUANTIFIER_ALL.to_string(),
        "filter".to_string(),
        OPERATOR_EQUAL_TO.to_string(),
        Value::String("PASS".to_string()),
        vec!["sample".to_string()]
    );
    variant_filters.push(variant_filter);

    variants_list.filter(variant_filters, 1);

    assert_eq!(variants_list.variants.len(), 1);
}

/// Verifies that `VariantsList::filter` function filters variants by `alternate_allele_read_count`.
///
/// # Scenario
/// Filter variants where alternate_allele_read_count < 3 for all variant calls.
///
/// # Expected
/// The variants list should have 1 variant.
#[test]
fn filter_variants_list_2() {
    let mut variants_list: VariantsList = VariantsList::new();
    let mut variant_1: Variant = Variant::new("variant_1".to_string());
    let mut variant_2: Variant = Variant::new("variant_2".to_string());
    let variant_call_1: VariantCall = VariantCall::new(
        "variant_call_1".to_string(),
        "sample".to_string(),
        "chr1".to_string(),
        100,
        "chr1".to_string(),
        100,
        INSERTION.to_string(),
        "T".to_string(),
        "TAAA".to_string(),
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
        3,
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
    variant_call_2.position_1 = 1000;
    variant_call_2.position_2 = 1000;
    variant_call_2.alternate_allele_read_count = 1;
    let mut variant_call_3: VariantCall = variant_call_1.clone();
    variant_call_3.id = "variant_call_3".to_string();
    variant_call_3.position_1 = 1000;
    variant_call_3.position_2 = 1000;
    variant_call_3.alternate_allele_read_count = 3;
    variant_call_3.variant_calling_method = "method2".to_string();
    variant_1.add_variant_call(variant_call_1);
    variant_2.add_variant_call(variant_call_2);
    variant_2.add_variant_call(variant_call_3);
    variants_list.add_variant(variant_1);
    variants_list.add_variant(variant_2);

    let mut variant_filters: Vec<VariantFilter> = Vec::new();
    let variant_filter: VariantFilter = VariantFilter::new(
        QUANTIFIER_ALL.to_string(),
        "alternate_allele_read_count".to_string(),
        OPERATOR_GREATER_THAN_EQUAL_TO.to_string(),
        Value::Number(3.into()),
        vec!["sample".to_string()]
    );
    variant_filters.push(variant_filter);

    variants_list.filter(variant_filters, 1);

    assert_eq!(variants_list.variants.len(), 1);
}

/// Verifies that `VariantsList::filter` function filters variants by `alternate_allele_read_count`.
///
/// # Scenario
/// Filter variants where alternate_allele_read_count < 3 for any variant call.
///
/// # Expected
/// The variants list should have 2 variants.
#[test]
fn filter_variants_list_3() {
    let mut variants_list: VariantsList = VariantsList::new();
    let mut variant_1: Variant = Variant::new("variant_1".to_string());
    let mut variant_2: Variant = Variant::new("variant_2".to_string());
    let variant_call_1: VariantCall = VariantCall::new(
        "variant_call_1".to_string(),
        "sample".to_string(),
        "chr1".to_string(),
        100,
        "chr1".to_string(),
        100,
        INSERTION.to_string(),
        "T".to_string(),
        "TAAA".to_string(),
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
        3,
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
    variant_call_2.position_1 = 1000;
    variant_call_2.position_2 = 1000;
    variant_call_2.alternate_allele_read_count = 1;
    let mut variant_call_3: VariantCall = variant_call_1.clone();
    variant_call_3.id = "variant_call_3".to_string();
    variant_call_3.position_1 = 1000;
    variant_call_3.position_2 = 1000;
    variant_call_3.alternate_allele_read_count = 3;
    variant_call_3.variant_calling_method = "method2".to_string();
    variant_1.add_variant_call(variant_call_1);
    variant_2.add_variant_call(variant_call_2);
    variant_2.add_variant_call(variant_call_3);
    variants_list.add_variant(variant_1);
    variants_list.add_variant(variant_2);

    let mut variant_filters: Vec<VariantFilter> = Vec::new();
    let variant_filter: VariantFilter = VariantFilter::new(
        QUANTIFIER_ANY.to_string(),
        "alternate_allele_read_count".to_string(),
        OPERATOR_GREATER_THAN_EQUAL_TO.to_string(),
        Value::Number(3.into()),
        vec!["sample".to_string()]
    );
    variant_filters.push(variant_filter);

    variants_list.filter(variant_filters, 1);

    assert_eq!(variants_list.variants.len(), 2);
}