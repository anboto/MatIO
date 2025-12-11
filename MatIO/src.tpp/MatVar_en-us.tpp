topic "MatVar";
[i448;a25;kKO9;2 $$1,0#37138531426314131252341829483380:class]
[l288;2 $$2,2#27521748481378242620020725143825:desc]
[0 $$3,0#96390100711032703541132217272105:end]
[H6;0 $$4,0#05600065144404261032431302351956:begin]
[i448;a25;kKO9;2 $$5,0#37138531426314131252341829483370:item]
[l288;a4;*@5;1 $$6,6#70004532496200323422659154056402:requirement]
[l288;i1121;b17;O9;~~~.1408;2 $$7,0#10431211400427159095818037425705:param]
[i448;b42;O9;2 $$8,8#61672508125594000341940100500538:tparam]
[b42;2 $$9,9#13035079074754324216151401829390:normal]
[2 $$0,0#00000000000000000000000000000000:Default]
[{_} 
[ {{10000@3 [s0;%% [*@(229)4 MatVar]]}}&]
[s2;^topic`:`/`/MatIO`/src`/MatVar`_en`-us^ &]
[s1;:Upp`:`:MatVar: [@(0.0.255) class]_[* MatVar]&]
[s2;%% &]
[s0;#%% MatVar is an RAII`-style wrapper around a matvar`_t`* that 
optionally owns the underlying resource and frees it on destruction. 
It provides convenient accessors for variable metadata (name, 
type, dimensions), total element counts, and struct field introspection 
and retrieval.&]
[s0;%% &]
[s3;%% &]
[ {{10000F(128)G(128)@1 [s0;%% [* Public Member List]]}}&]
[s4; &]
[s5;:Upp`:`:MatVar`:`:operator matvar`_t`*`(`): operator matvar`_t 
[@(0.0.255) `*]()&]
[s2;%% Returns the underlying matvar`_t`*, allowing interoperation 
with MatIO functions.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatVar`:`:GetName`(`): [@(0.0.255) const] [@(0.0.255) char] 
[@(0.0.255) `*][* GetName]()&]
[s2;%% Returns the variable’s name; this method throws if the variable 
has not been created.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatVar`:`:GetType`(`): [@(0.0.255) enum] matio`_classes 
[* GetType]()&]
[s2;%% Returns the MatIO class type (such as MAT`_C`_DOUBLE or MAT`_C`_STRUCT); 
this method throws if the variable is not valid.&]
[s3;2%% &]
[s4; &]
[s5;:Upp`:`:MatVar`:`:GetDimCount`(`)const: [@(0.0.255) int] [* GetDimCount]() 
[@(0.0.255) const]&]
[s2;%% Returns the number of dimensions (rank) of the variable; this 
method throws if the variable is not valid.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatVar`:`:GetDimCount`(int`)const: [@(0.0.255) int] [* GetDimCount]([@(0.0.255) i
nt] [*@3 dim]) [@(0.0.255) const]&]
[s2; [@N Returns the size of the specified dimension; this method throws 
if the variable is not valid.]&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatVar`:`:GetTypeString`(`): [@(0.0.255) const] [@(0.0.255) char] 
[@(0.0.255) `*][* GetTypeString]()&]
[s2;%% Returns a human`-readable description of the MatIO class type; 
this method throws if the variable is not valid.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatVar`:`:GetCount`(`): [@(0.0.255) int] [* GetCount]()&]
[s2;%% Returns the total number of elements, computed as the product 
of all dimensions; this method throws if the variable is not 
valid.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatVar`:`:GetFieldCount`(`): [@(0.0.255) int] [* GetFieldCount]()&]
[s2;%% Returns the number of fields in a struct; this method throws 
if the variable is not valid or not a struct.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatVar`:`:GetFieldName`(int`): [@(0.0.255) const] [@(0.0.255) char] 
[@(0.0.255) `*][* GetFieldName]([@(0.0.255) int] [*@3 id])&]
[s2; [@N Returns the name of a struct field by index; this method throws 
if the variable is not valid or not a struct.]&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatVar`:`:IsLoaded`(`): [@(0.0.255) bool] [* IsLoaded]()&]
[s2;%% Returns true if the wrapper currently holds a non`-null variable 
pointer.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatVar`:`:IsStruct`(`): [@(0.0.255) bool] [* IsStruct]()&]
[s2;%% Returns true if the variable’s class type indicates a MATLAB 
struct.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatVar`:`:GetStructVar`(const char`*`): MatVar [* GetStructVar]([@(0.0.255) c
onst] [@(0.0.255) char] [@(0.0.255) `*][*@3 field`_name])&]
[s2; [@N Returns a wrapper for the named struct field without taking 
ownership of the returned matvar`_t`*.]&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatVar`:`:ExistStructVar`(const char`*`): [@(0.0.255) bool] 
[* ExistStructVar]([@(0.0.255) const] [@(0.0.255) char] [@(0.0.255) `*][*@3 field`_name])&]
[s2;%% Returns true if the variable [%-*@3 field`_name ]exists.&]
[s3; &]
[ {{10000F(128)G(128)@1 [s0;%% [* Constructor Detail]]}}&]
[s4; &]
[s5;:Upp`:`:MatVar`:`:MatVar`(`): [* MatVar]()&]
[s2;%% Constructs an empty wrapper that does not reference any variable&]
[s3;%% &]
[s4; &]
[s5;:Upp`:`:MatVar`:`:MatVar`(matvar`_t`*`,bool`): [* MatVar](matvar`_t 
[@(0.0.255) `*][*@3 var`_], [@(0.0.255) bool] [*@3 del`_] [@(0.0.255) `=] 
[@(0.0.255) true])&]
[s2;%% Wraps an existing variable pointer and optionally assumes 
ownership for freeing it.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatVar`:`:`~`(`): [* `~MatVar]()&]
[s2;%% Frees the underlying matvar`_t if the wrapper owns it .&]
[s3; &]
[s0; ]]