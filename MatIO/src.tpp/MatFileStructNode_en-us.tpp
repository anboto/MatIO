topic "";
[l288;2 $$1,1#27521748481378242620020725143825:desc]
[i448;a25;kKO9;2 $$2,0#37138531426314131252341829483380:class]
[0 $$3,0#96390100711032703541132217272105:end]
[H6;0 $$4,0#05600065144404261032431302351956:begin]
[i448;a25;kKO9;2 $$5,0#37138531426314131252341829483370:item]
[ $$0,0#00000000000000000000000000000000:Default]
[{_} 
[ {{10000@3 [s0;%% [*@(229)4 MatFile`::StructNode]]}}&]
[s1; &]
[s2;:Upp`:`:MatFile`:`:StructNode: [@(0.0.255) class]_[* MatFile`::StructNode]&]
[s1;%% &]
[s0;#%% StructNode is a builder for creating and writing a 1×1 MATLAB 
struct backed by storage owned by the node itself. It keeps all 
field data alive until the node is written to the associated 
MatFile, which ensures that the MatIO library receives valid 
pointers during serialization. It supports scalar fields, string 
fields, numeric and complex matrices, and nested child structs.&]
[s0;%% &]
[s3;%% &]
[ {{10000F(128)G(128)@1 [s0;%% [* Public Member List]]}}&]
[s4; &]
[s5;:Upp`:`:MatFile`:`:StructNode`:`:StructNode`(`): [* StructNode]()&]
[s1;%% Default constructor.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:StructNode`:`:StructNode`(MatFile`&`,const char`*`,const Vector`&`): [* S
tructNode](MatFile[@(0.0.255) `&] [*@3 mf], [@(0.0.255) const] [@(0.0.255) char] 
[@(0.0.255) `*][*@3 name], [@(0.0.255) const] Vector<String>[@(0.0.255) `&] 
[*@3 fieldNames])&]
[s1; [@N Constructor with file association and predefined fields.]&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:StructNode`:`:StructNode`(const StructNode`&`)`=delete: [* Struc
tNode]([@(0.0.255) const] [* StructNode][@(0.0.255) `&])&]
[s1;%% Copy operations are disabled to avoid accidental sharing of 
owned matvar`_t`*.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:StructNode`:`:StructNode`(StructNode`&`&`): [* StructNode]([* Stru
ctNode][@(0.0.255) `&`&] other)&]
[s1;%% Move operations transfer ownership of the underlying matvar`_t`* 
and backing storage.&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:StructNode`:`:`~`(`): [* `~StructNode]()&]
[s1;%% Frees the internal matvar`_t`* unless it has been transferred 
to the file by Write.&]
[s3; &]
[s3; &]
[ {{10000F(128)G(128)@1 [s0;%% [* Constructor Detail]]}}&]
[s3; &]
[s5;:Upp`:`:MatFile`:`:StructNode`:`:Create`(const char`*`,const Vector`&`): [@(0.0.255) v
oid] [* Create]([@(0.0.255) const] [@(0.0.255) char] [@(0.0.255) `*][*@3 name], 
[@(0.0.255) const] Vector<String>[@(0.0.255) `&] [*@3 fieldNames])&]
[s1; [@N Creates or recreates a struct with the specified name and 
fields.]&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:StructNode`:`:Write`(const char`*`,int`): [@(0.0.255) void] 
[* Write]([@(0.0.255) const] [@(0.0.255) char] [@(0.0.255) `*][*@3 field], 
[@(0.0.255) int] [*@3 v])&]
[s1; [@N Adds an integer scalar field to the struct.]&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:StructNode`:`:Write`(const char`*`,double`): [@(0.0.255) void] 
[* Write]([@(0.0.255) const] [@(0.0.255) char] [@(0.0.255) `*][*@3 field], 
[@(0.0.255) double] [*@3 v])&]
[s1; [@N Adds a double scalar field to the struct.]&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:StructNode`:`:Write`(const char`*`,const String`&`): [@(0.0.255) v
oid] [* Write]([@(0.0.255) const] [@(0.0.255) char] [@(0.0.255) `*][*@3 field], 
[@(0.0.255) const] String[@(0.0.255) `&] [*@3 s])&]
[s1; [@N Adds a character vector field, stored as UTF‑8 when available 
or uint8 otherwise.]&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:StructNode`:`:WriteColMajor`(const char`*`,size`_t`,size`_t`,const double`*`): [@(0.0.255) v
oid] [* WriteColMajor]([@(0.0.255) const] [@(0.0.255) char] [@(0.0.255) `*][*@3 field], 
size`_t [*@3 rows], size`_t [*@3 cols], [@(0.0.255) const] [@(0.0.255) double] 
[@(0.0.255) `*][*@3 colMajor])&]
[s1; [@N Adds a double matrix field (column`-major input).]&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:StructNode`:`:WriteColMajor`(const char`*`,size`_t`,size`_t`,const int`*`): [@(0.0.255) v
oid] [* WriteColMajor]([@(0.0.255) const] [@(0.0.255) char] [@(0.0.255) `*][*@3 field], 
size`_t [*@3 rows], size`_t [*@3 cols], [@(0.0.255) const] [@(0.0.255) int] 
[@(0.0.255) `*][*@3 colMajor])&]
[s1; [@N Adds an integer matrix field (column`-major input).]&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:StructNode`:`:AddChild`(const char`*`,const Vector`&`): Struct
Node[@(0.0.255) `&] [* AddChild]([@(0.0.255) const] [@(0.0.255) char] 
[@(0.0.255) `*][*@3 fieldName], [@(0.0.255) const] Vector<String>[@(0.0.255) `&] 
[*@3 childFields])&]
[s1; [@N Creates a child struct as a field of the current struct and 
returns a reference to the child for population.]&]
[s3; &]
[s4; &]
[s5;:Upp`:`:MatFile`:`:StructNode`:`:Write`(matio`_compression`): [@(0.0.255) void] 
[* Write](matio`_compression compression [@(0.0.255) `=] [*@3 MAT`_COMPRESSION`_NONE])&]
[s1; [@N Finalizes the struct, writes it to the associated MatFile, 
and relinquishes ownership to MatIO.]&]
[s3; &]
[s0; ]]