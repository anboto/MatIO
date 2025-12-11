topic "MatFile";
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
[ {{10000@3 [s0;%% [*@(229)4 MatFile]]}}&]
[s2; &]
[s1;:Upp`:`:MatFile: [@(0.0.255) class]_[* MatFile]&]
[s2;%% &]
[s0;#%% MatFile is an RAII wrapper that opens, creates, reads, and 
writes MATLAB .mat files using the MatIO library. It supports 
MAT4, MAT5, and MAT7.3 formats, and it provides convenience methods 
for writing scalar values, strings, matrices (double, int, and 
complex), U`+`+ containers, and Eigen matrices. It also includes 
static reading helpers for extracting data in row`-major buffers 
and Eigen`-friendly forms, plus a nested StructNode for building 
and writing structured variables.&]
[s0;%% &]
[s3;%% &]
[ {{10000F(128)G(128)@1 [s0;%% [* Public Member List]]}}&]
[s4; &]
[s5;:Upp`:`:MatFile`:`:OpenRead`(const char`*`): [@(0.0.255) bool] 
[* OpenRead]([@(0.0.255) const] [@(0.0.255) char] [@(0.0.255) `*][*@3 path])&]
[s2;%% Opens a .mat file read`-only; returns true if the file handle 
is valid.&]
[s3;%% &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:OpenReadWrite`(const char`*`): [@(0.0.255) bool] 
[* OpenReadWrite]([@(0.0.255) const] [@(0.0.255) char] [@(0.0.255) `*][*@3 path])&]
[s2;%% Opens a .mat file for read and write; returns true if the 
file handle is valid.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:OpenCreate`(const char`*`,mat`_ft`): [@(0.0.255) bool] 
[* OpenCreate]([@(0.0.255) const] [@(0.0.255) char] [@(0.0.255) `*][*@3 path], 
mat`_ft ver [@(0.0.255) `=] [@(0.0.255) enum] mat`_ft[@(0.0.255) `::][*@3 MAT`_FT`_MAT5])&]
[s2; [@N Creates a new .mat file for the specified format version; 
returns true if the file handle is valid.]&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:GetVersion`(`): mat`_ft [* GetVersion]()&]
[s2;%% Returns the format version of the open file; this method throws 
if the file is not open.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:GetVersionName`(`): String [* GetVersionName]()&]
[s2;%% Returns a human`-readable version string (`"4`", `"5`", `"7.3`", 
or `"unknown`").&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:GetVarList`(`)const: Vector<String> [* GetVarList]() 
[@(0.0.255) const]&]
[s2;%% Returns the list of top`-level variable names in the open 
file.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:GetVarName`(const char`*`)const: String [* GetVarName]([@(0.0.255) c
onst] [@(0.0.255) char] [@(0.0.255) `*][*@3 name]) [@(0.0.255) const]&]
[s2;%% Returns the real name of a variable, independently of its 
case.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:Mat`_t`(`)const: mat`_t [@(0.0.255) `*][* Mat`_t]() 
[@(0.0.255) const]&]
[s2;%% Returns the raw mat`_t`* file handle for low`-level operations.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:GetVar`(const char`*`,bool`)const: MatVar 
[* GetVar]([@(0.0.255) const] [@(0.0.255) char] [@(0.0.255) `*][*@3 name], 
[@(0.0.255) bool] [*@3 nocase] [@(0.0.255) `=] [@(0.0.255) false]) [@(0.0.255) const]&]
[s2;%% [%-@N Reads a top`-level variable by name and returns an owning 
MatVar; this method throws if the file is not open or the variable 
cannot be found. ][%-*@3 nocase][%-  ][%-@N allows to input variable 
name ]independently of its case.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:GetVar`(const Vector`&`,bool`)const: MatVar 
[* GetVar]([@(0.0.255) const] Vector<String>[@(0.0.255) `&] [*@3 names], 
[@(0.0.255) bool] [*@3 nocase] [@(0.0.255) `=] [@(0.0.255) false]) [@(0.0.255) const]&]
[s2;%% [%-@N Reads a top`-level variable from a list of synonym names 
and returns an owning MatVar; this method throws if the file 
is not open or the variable cannot be found. ][%-*@3 nocase][%-  
][%-@N allows to input variable name ]independently of its case.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:Exist`(const char`*`)const: [@(0.0.255) bool] 
[* Exist]([@(0.0.255) const] [@(0.0.255) char] [@(0.0.255) `*][*@3 name]) 
[@(0.0.255) const]&]
[s2;%% Returns true if the variable  [%-*@3 name ]exists.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:Write`(const char`*`,int`,matio`_compression`): [@(0.0.255) void
] [* Write]([@(0.0.255) const] [@(0.0.255) char] [@(0.0.255) `*][*@3 name], 
[@(0.0.255) int] [*@3 v], matio`_compression compression [@(0.0.255) `=] 
[*@3 MAT`_COMPRESSION`_NONE])&]
[s2; [@N Writes an integer scalar to the file.]&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:Write`(const char`*`,double`,matio`_compression`): [@(0.0.255) v
oid] [* Write]([@(0.0.255) const] [@(0.0.255) char] [@(0.0.255) `*][*@3 name], 
[@(0.0.255) double] [*@3 v], matio`_compression compression [@(0.0.255) `=] 
[*@3 MAT`_COMPRESSION`_NONE])&]
[s2; [@N Writes a double`-precision scalar to the file.]&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:Write`(const char`*`,const char`*`,matio`_compression`): [@(0.0.255) v
oid] [* Write]([@(0.0.255) const] [@(0.0.255) char] [@(0.0.255) `*][*@3 name], 
[@(0.0.255) const] [@(0.0.255) char] [@(0.0.255) `*][*@3 s], matio`_compression 
compression [@(0.0.255) `=] [*@3 MAT`_COMPRESSION`_NONE])&]
[s2; [@N Writes a character vector using UTF‑8 (when available) or 
as uint8.]&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:WriteColMajor`(const char`*`,size`_t`,size`_t`,const double`*`,matio`_compression`): [@(0.0.255) v
oid] [* WriteColMajor]([@(0.0.255) const] [@(0.0.255) char] [@(0.0.255) `*][*@3 name], 
size`_t [*@3 rows], size`_t [*@3 cols], [@(0.0.255) const] [@(0.0.255) double] 
[@(0.0.255) `*][*@3 colMajor], matio`_compression compression [@(0.0.255) `=] 
[*@3 MAT`_COMPRESSION`_NONE])&]
[s2; [@N Writes a double matrix from column`-major memory.]&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:WriteColMajor`(const char`*`,size`_t`,size`_t`,const int`*`,matio`_compression`): [@(0.0.255) v
oid] [* WriteColMajor]([@(0.0.255) const] [@(0.0.255) char] [@(0.0.255) `*][*@3 name], 
size`_t [*@3 rows], size`_t [*@3 cols], [@(0.0.255) const] [@(0.0.255) int] 
[@(0.0.255) `*][*@3 colMajor], matio`_compression compression [@(0.0.255) `=] 
[*@3 MAT`_COMPRESSION`_NONE])&]
[s2;%%  [%-*@3 name] [%-*@3 rows] [%-*@3 cols] [%-*@3 colMajor] [%-*@3 MAT`_COMPRESSION`_NONE] 
.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:WriteRowMajor`(const char`*`,size`_t`,size`_t`,const T`*`,matio`_compression`): [@(0.0.255) t
emplate] <[@(0.0.255) class] T> [@(0.0.255) void] [* WriteRowMajor]([@(0.0.255) const] 
[@(0.0.255) char] [@(0.0.255) `*][*@3 name], size`_t [*@3 rows], size`_t 
[*@3 cols], [@(0.0.255) const] T [@(0.0.255) `*][*@3 rowMajor], matio`_compression 
compression [@(0.0.255) `=] [*@3 MAT`_COMPRESSION`_NONE])&]
[s2; [@N  Writes an integer matrix from column`-major memory.]&]
[s3; &]
[s3;%% &]
[ {{10000F(128)G(128)@1 [s0;%% [* Constructor Detail]]}}&]
[s3; &]
[s5;:Upp`:`:MatFile`:`:MatFile`(`): [* MatFile]()&]
[s2; Constructs an empty file wrapper.&]
[s3; &]
[s4;%% &]
[s5;:Upp`:`:MatFile`:`:`~`(`): [* `~MatFile]()&]
[s2;%% Closes the file handle if it is open.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:MatFile`(const MatFile`&`)`=delete: [* MatFile]([@(0.0.255) const] 
[* MatFile][@(0.0.255) `&])&]
[s2;%% Copy operations are disabled to avoid double`-close or shared 
handle issues.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:MatFile`(MatFile`&`&`): [* MatFile]([* MatFile][@(0.0.255) `&`&] 
o)&]
[s2;%% Move operations transfer ownership of the file handle and 
path..&]
[s3; &]
[s0;%% ]]